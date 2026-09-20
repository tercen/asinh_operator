//! asinh_operator — Rust port of `tercen/asinh_operator` (fixed and manual methods).
//!
//! The crosstab contract: rows are channels, columns are events, `y` is the value, and the
//! result is `asinh(y / cofactor)` per cell, carried back on `.ri` / `.ci`. `fixed` takes the
//! cofactor from the `scale` property; `manual` takes one per channel from the second row
//! factor. `auto` (flowVS estimation) is deliberately not here — see `props::settings_from_ctx`.
pub mod algorithm;
#[cfg(feature = "auto")]
pub mod cofactors;
pub mod input;
pub mod output;
pub mod pagecache;
pub mod progress;
pub mod props;
pub mod tson;
pub mod upload;

use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use anyhow::{Context, Result, anyhow};
use tercen_rs::context::ContextBase;
use tercen_rs::{DevContext, ProductionContext, TercenClient};

use progress::Reporter;
use props::{Method, Settings};
use tson::TsonWriter;

/// Cells fetched per gRPC round trip.
const CHUNK: i64 = 1_000_000;

pub fn init_tracing() {
    let filter = tracing_subscriber::EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info"));
    let _ = tracing_subscriber::fmt().with_env_filter(filter).try_init();
}

pub fn require_env(name: &str) -> Result<String> {
    std::env::var(name).map_err(|_| anyhow!("{name} is not set"))
}

/// Production entry point (`--taskId`).
pub async fn run(task_id: &str) -> Result<()> {
    tracing::info!("asinh_operator starting (task_id={task_id})");
    let client = build_client().await?;
    let ctx = ProductionContext::from_task_id(client, task_id)
        .await
        .map_err(|e| anyhow!("load task {task_id}: {e}"))?;
    execute(
        &ctx,
        Mode::Production {
            task_id: task_id.to_string(),
        },
    )
    .await
}

/// Dev entry point (`WORKFLOW_ID` / `STEP_ID`).
pub async fn run_dev(workflow_id: &str, step_id: &str) -> Result<()> {
    tracing::info!("asinh_operator starting in dev mode ({workflow_id} / {step_id})");
    let client = build_client().await?;
    let ctx = DevContext::from_workflow_step(client, workflow_id, step_id)
        .await
        .map_err(|e| anyhow!("load workflow {workflow_id} / step {step_id}: {e}"))?;
    execute(
        &ctx,
        Mode::Dev {
            workflow_id: workflow_id.to_string(),
            step_id: step_id.to_string(),
        },
    )
    .await
}

enum Mode {
    Production {
        task_id: String,
    },
    Dev {
        workflow_id: String,
        step_id: String,
    },
}

async fn build_client() -> Result<Arc<TercenClient>> {
    let client = TercenClient::from_env()
        .await
        .map_err(|e| anyhow!("connect to Tercen: {e}"))?;
    tracing::info!("connected to Tercen");
    Ok(Arc::new(client))
}

async fn execute(ctx: &ContextBase, mode: Mode) -> Result<()> {
    let t_start = Instant::now();
    tracing::info!(
        workflow = ctx.workflow_id(),
        step = ctx.step_id(),
        namespace = ctx.namespace(),
        "context loaded"
    );
    let rep = match &mode {
        Mode::Production { task_id } => Reporter::spawn(Arc::clone(ctx.client()), task_id.clone()),
        Mode::Dev { .. } => Reporter::silent(),
    };
    let s = props::settings_from_ctx(ctx)?;
    tracing::info!(?s, "properties");

    #[allow(unused_mut)] // only the `auto` path fills it
    let mut cofactor_table: Vec<output::CofactorRow> = Vec::new();
    let cofactors = match s.method {
        Method::Fixed => None,
        Method::Manual => {
            let v = input::row_cofactors(ctx).await?;
            tracing::info!(channels = v.len(), "per-channel cofactors read");
            Some(v)
        }
        #[cfg(feature = "auto")]
        Method::Auto => {
            rep.at(0, "Estimating cofactors (flowVS)");
            let est = cofactors::estimate(ctx, &s).await?;
            for row in &est.table {
                rep.info(format!(
                    "{}: cofactor {:.1}{}",
                    row.channel,
                    row.cofactor,
                    if row.unstable { " (unstable)" } else { "" }
                ));
            }
            cofactor_table = est.table;
            Some(est.per_row)
        }
    };
    let with_cofactors = !cofactor_table.is_empty();

    let n_cells = input::cell_count(ctx).await?;
    let collect = n_cells <= output::COLLECT_MAX_CELLS;
    tracing::info!(
        n_cells,
        mode = if collect { "collect" } else { "stream" },
        "crosstab size"
    );
    rep.at(0, format!("Transforming {n_cells} values"));

    let work_root =
        std::env::temp_dir().join(format!("asinh_op_{}_{}", ctx.workflow_id(), ctx.step_id()));
    std::fs::create_dir_all(&work_root)
        .with_context(|| format!("create {}", work_root.display()))?;
    let _guard = TempDirGuard(work_root.clone());
    let result_path = work_root.join("result.tson");

    let t = Instant::now();
    {
        let f = std::fs::File::create(&result_path)
            .with_context(|| format!("create {}", result_path.display()))?;
        // The result is one row per cell and can be GBs; hand written pages back to the kernel
        // or the cgroup counts them against --memory (pagecache.rs).
        let w = std::io::BufWriter::with_capacity(4 << 20, pagecache::Releasing::new(f, 256 << 20));
        let mut w = TsonWriter::new(w)?;
        if collect {
            write_collect(
                ctx,
                &s,
                cofactors.as_deref(),
                n_cells,
                &mut w,
                &rep,
                with_cofactors,
            )
            .await?;
        } else {
            write_streamed(
                ctx,
                &s,
                cofactors.as_deref(),
                n_cells,
                &mut w,
                &rep,
                with_cofactors,
            )
            .await?;
        }
        if with_cofactors {
            output::write_cofactor_table(&mut w, &cofactor_table)?;
        }
        output::write_footer(&mut w, with_cofactors)?;
    }
    let bytes = std::fs::metadata(&result_path)?.len();
    let secs = t.elapsed().as_secs_f64();
    tracing::info!(
        bytes,
        secs = format!("{secs:.1}"),
        cells_per_s = format!("{:.0}", n_cells as f64 / secs.max(1e-9)),
        "result written"
    );
    pagecache::release_path(&result_path);

    rep.at(progress::UPLOAD.0, "Uploading the result");
    match mode {
        Mode::Production { task_id } => {
            upload::save_production(ctx, &task_id, &result_path, &rep).await?
        }
        Mode::Dev {
            workflow_id,
            step_id,
        } => {
            let saved = upload::save_dev(ctx, &workflow_id, &step_id, &result_path).await?;
            tracing::info!(
                task_id = saved.task_id,
                file_id = saved.file_id,
                "dev result saved"
            );
        }
    }
    rep.at(100, "Done");
    rep.info(format!(
        "asinh complete: {n_cells} values in {:.1} s",
        t_start.elapsed().as_secs_f64()
    ));
    tracing::info!(
        total_secs = format!("{:.1}", t_start.elapsed().as_secs_f64()),
        peak_rss_kb = peak_rss_kb().unwrap_or(0),
        "done"
    );
    Ok(())
}

/// One pass, buffered. Used under `COLLECT_MAX_CELLS`.
async fn write_collect<W: std::io::Write>(
    ctx: &ContextBase,
    s: &Settings,
    cofactors: Option<&[f64]>,
    n_cells: usize,
    w: &mut TsonWriter<W>,
    rep: &Reporter,
    with_cofactors: bool,
) -> Result<()> {
    let mut ri: Vec<i32> = Vec::with_capacity(n_cells);
    let mut ci: Vec<i32> = Vec::with_capacity(n_cells);
    let mut y: Vec<f64> = Vec::with_capacity(n_cells);
    let mut seen = 0usize;
    let mut err: Option<anyhow::Error> = None;
    ctx.streamer()
        .stream_table_chunked(
            ctx.qt_hash(),
            Some(vec![".ri".into(), ".ci".into(), ".y".into()]),
            CHUNK,
            |bytes| {
                match input::decode_chunk(&bytes, true, true, true) {
                    Ok(c) => {
                        seen += c.y.len();
                        ri.extend_from_slice(&c.ri);
                        ci.extend_from_slice(&c.ci);
                        y.extend_from_slice(&c.y);
                        rep.at(
                            progress::band(progress::READ, seen, n_cells.max(1)),
                            format!("Read {seen} of {n_cells} values"),
                        );
                    }
                    Err(e) => err = Some(e),
                }
                Ok(())
            },
        )
        .await
        .map_err(|e| anyhow!("stream the crosstab: {e}"))?;
    if let Some(e) = err {
        return Err(e);
    }
    transform(&mut y, &ri, s, cofactors)?;

    let ns = output::value_column(ctx.namespace());
    let cols = output::result_columns(&ns);
    let n = y.len();
    output::write_header(
        w,
        &uuid_like(ctx),
        n,
        &cols,
        1 + usize::from(with_cofactors),
    )?;
    output::write_column_header(w, &cols[0], n)?;
    w.i32_list(&ri)?;
    output::write_column_header(w, &cols[1], n)?;
    w.i32_list(&ci)?;
    output::write_column_header(w, &cols[2], n)?;
    rep.at(progress::WRITE.0, "Writing the result");
    w.f64_list(&y)?;
    Ok(())
}

/// Three passes, one per column, constant memory. Used above `COLLECT_MAX_CELLS`.
async fn write_streamed<W: std::io::Write>(
    ctx: &ContextBase,
    s: &Settings,
    cofactors: Option<&[f64]>,
    n_cells: usize,
    w: &mut TsonWriter<W>,
    rep: &Reporter,
    with_cofactors: bool,
) -> Result<()> {
    let ns = output::value_column(ctx.namespace());
    let cols = output::result_columns(&ns);
    output::write_header(
        w,
        &uuid_like(ctx),
        n_cells,
        &cols,
        1 + usize::from(with_cofactors),
    )?;

    for (pass, col) in [".ri", ".ci"].iter().enumerate() {
        output::write_column_header(w, &cols[pass], n_cells)?;
        w.i32_list_header(n_cells)?;
        let mut written = 0usize;
        let mut err: Option<anyhow::Error> = None;
        ctx.streamer()
            .stream_table_chunked(ctx.qt_hash(), Some(vec![col.to_string()]), CHUNK, |bytes| {
                match input::decode_chunk(&bytes, *col == ".ri", *col == ".ci", false) {
                    Ok(c) => {
                        let v = if *col == ".ri" { &c.ri } else { &c.ci };
                        if let Err(e) = w.i32_chunk(v) {
                            err = Some(e.into());
                        }
                        written += v.len();
                        rep.at(
                            progress::band(progress::READ, written, n_cells.max(1)),
                            format!("Pass {}: {written} of {n_cells}", pass + 1),
                        );
                    }
                    Err(e) => err = Some(e),
                }
                Ok(())
            })
            .await
            .map_err(|e| anyhow!("stream {col}: {e}"))?;
        if let Some(e) = err {
            return Err(e);
        }
        check_count(written, n_cells, col)?;
    }

    output::write_column_header(w, &cols[2], n_cells)?;
    w.f64_list_header(n_cells)?;
    let mut written = 0usize;
    let mut err: Option<anyhow::Error> = None;
    ctx.streamer()
        .stream_table_chunked(
            ctx.qt_hash(),
            Some(vec![".ri".into(), ".y".into()]),
            CHUNK,
            |bytes| {
                match input::decode_chunk(&bytes, true, false, true) {
                    Ok(mut c) => {
                        if let Err(e) = transform(&mut c.y, &c.ri, s, cofactors) {
                            err = Some(e);
                            return Ok(());
                        }
                        if let Err(e) = w.f64_chunk(&c.y) {
                            err = Some(e.into());
                        }
                        written += c.y.len();
                        rep.at(
                            progress::band(progress::WRITE, written, n_cells.max(1)),
                            format!("Transformed {written} of {n_cells} values"),
                        );
                    }
                    Err(e) => err = Some(e),
                }
                Ok(())
            },
        )
        .await
        .map_err(|e| anyhow!("stream values: {e}"))?;
    if let Some(e) = err {
        return Err(e);
    }
    check_count(written, n_cells, ".y")?;
    Ok(())
}

fn transform(y: &mut [f64], ri: &[i32], s: &Settings, cofactors: Option<&[f64]>) -> Result<()> {
    match cofactors {
        None => algorithm::asinh_fixed(y, s.scale),
        Some(c) => algorithm::asinh_per_row(y, ri, c).map_err(|i| {
            anyhow!(
                "row index {i} has no cofactor — the cofactor row factor has {} values, which \
                 does not cover the projected rows",
                c.len()
            )
        })?,
    }
    Ok(())
}

/// The streamed writer declares `nRows` before it has the data, so a short or long stream would
/// silently corrupt the result. Fail instead.
fn check_count(written: usize, expected: usize, what: &str) -> Result<()> {
    if written != expected {
        return Err(anyhow!(
            "the crosstab returned {written} values for {what} but its schema says {expected}; \
             the table changed under the operator"
        ));
    }
    Ok(())
}

/// A stable per-run table name (`save_table` uses a uuid; the value is never read back).
fn uuid_like(ctx: &ContextBase) -> String {
    format!("{}_{}", ctx.step_id(), ctx.qt_hash())
}

fn peak_rss_kb() -> Option<u64> {
    let s = std::fs::read_to_string("/proc/self/status").ok()?;
    s.lines()
        .find(|l| l.starts_with("VmHWM:"))?
        .split_whitespace()
        .nth(1)?
        .parse()
        .ok()
}

struct TempDirGuard(PathBuf);
impl Drop for TempDirGuard {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
