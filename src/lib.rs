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
const CHUNK: usize = 1_000_000;

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

    let n_cells = input::cell_count(ctx).await?;
    let collect = n_cells <= output::COLLECT_MAX_CELLS;
    #[allow(unused_mut)] // only the `auto` path fills it
    let mut cofactor_table: Vec<output::CofactorRow> = Vec::new();
    #[cfg(feature = "auto")]
    let with_cofactors = s.method == Method::Auto;
    #[cfg(not(feature = "auto"))]
    let with_cofactors = false;

    let cofactors = match s.method {
        Method::Fixed => None,
        Method::Manual => {
            let v = input::row_cofactors(ctx).await?;
            tracing::info!(channels = v.len(), "per-channel cofactors read");
            Some(v)
        }
        // Collect mode estimates from the values it is about to buffer, so only the streaming
        // path spends a pass on it.
        #[cfg(feature = "auto")]
        Method::Auto if collect => None,
        #[cfg(feature = "auto")]
        Method::Auto => {
            rep.at(0, "Estimating cofactors (flowVS)");
            let est = cofactors::estimate(ctx, &s).await?;
            report_cofactors(&rep, &est.table);
            cofactor_table = est.table;
            Some(est.per_row)
        }
    };
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
            let table = write_collect(
                ctx,
                &s,
                cofactors.as_deref(),
                n_cells,
                &mut w,
                &rep,
                with_cofactors,
            )
            .await?;
            if !table.is_empty() {
                report_cofactors(&rep, &table);
                cofactor_table = table;
            }
        } else {
            write_streamed(
                ctx,
                &s,
                cofactors.as_deref(),
                n_cells,
                &mut w,
                &rep,
                with_cofactors,
                &work_root,
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
) -> Result<Vec<output::CofactorRow>> {
    let mut ri: Vec<i32> = Vec::with_capacity(n_cells);
    let mut ci: Vec<i32> = Vec::with_capacity(n_cells);
    let mut y: Vec<f64> = Vec::with_capacity(n_cells);
    let mut seen = 0usize;
    input::for_each_chunk(ctx, &[".ri", ".ci", ".y"], n_cells, CHUNK, |c| {
        seen += c.len();
        ri.extend_from_slice(&c.ri);
        ci.extend_from_slice(&c.ci);
        y.extend_from_slice(&c.y);
        rep.at(
            progress::band(progress::READ, seen, n_cells.max(1)),
            format!("Read {seen} of {n_cells} values"),
        );
        Ok(())
    })
    .await?;

    // `auto` costs no extra pass here: the values are already in memory, so the subsample comes
    // from them instead of reading the crosstab a second time.
    #[cfg(feature = "auto")]
    let estimated = if s.method == Method::Auto {
        Some(cofactors::estimate_from_cells(ctx, s, &ri, &ci, &y).await?)
    } else {
        None
    };
    #[cfg(feature = "auto")]
    let owned: Option<Vec<f64>> = estimated.as_ref().map(|e| e.per_row.clone());
    #[cfg(feature = "auto")]
    let cofactors = owned.as_deref().or(cofactors);
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
    rep.at(progress::WRITE.0, "Writing the result");
    output::write_column_header(w, &cols[0], n)?;
    w.f64_list(&y)?;
    output::write_column_header(w, &cols[1], n)?;
    w.i32_list(&ri)?;
    output::write_column_header(w, &cols[2], n)?;
    w.i32_list(&ci)?;
    #[cfg(feature = "auto")]
    return Ok(estimated.map(|e| e.table).unwrap_or_default());
    #[cfg(not(feature = "auto"))]
    Ok(Vec::new())
}

/// One pass over the crosstab, with the two index columns spilled to disk.
///
/// TSON is column-major, so a column must be finished before the next begins, and the naive way
/// to do that without buffering is one pass per column. Transfer dominates this operator by
/// roughly twenty to one against the arithmetic, so three passes would be the most expensive
/// thing it does. Instead the value column is declared **first** and streamed straight through
/// while `.ri` and `.ci` go to temporary files in exactly the little-endian layout TSON wants;
/// the two files are then poured into the result. One pass, one chunk of memory, and 8 bytes per
/// cell of scratch disk that is handed back to the kernel as it goes.
#[allow(clippy::too_many_arguments)]
async fn write_streamed<W: std::io::Write>(
    ctx: &ContextBase,
    s: &Settings,
    cofactors: Option<&[f64]>,
    n_cells: usize,
    w: &mut TsonWriter<W>,
    rep: &Reporter,
    with_cofactors: bool,
    work_root: &std::path::Path,
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

    let ri_path = work_root.join("ri.i32");
    let ci_path = work_root.join("ci.i32");
    let mut written = 0usize;
    {
        let mut ri_out = std::io::BufWriter::with_capacity(
            1 << 20,
            std::fs::File::create(&ri_path)
                .with_context(|| format!("create {}", ri_path.display()))?,
        );
        let mut ci_out = std::io::BufWriter::with_capacity(
            1 << 20,
            std::fs::File::create(&ci_path)
                .with_context(|| format!("create {}", ci_path.display()))?,
        );
        output::write_column_header(w, &cols[0], n_cells)?;
        w.f64_list_header(n_cells)?;

        input::for_each_chunk(ctx, &[".ri", ".ci", ".y"], n_cells, CHUNK, |mut c| {
            transform(&mut c.y, &c.ri, s, cofactors)?;
            w.f64_chunk(&c.y)?;
            write_i32_le(&mut ri_out, &c.ri)?;
            write_i32_le(&mut ci_out, &c.ci)?;
            written += c.len();
            rep.at(
                progress::band(progress::READ, written, n_cells.max(1)),
                format!("Transformed {written} of {n_cells} values"),
            );
            Ok(())
        })
        .await?;
        use std::io::Write as _;
        ri_out.flush()?;
        ci_out.flush()?;
    }
    check_count(written, n_cells, "the crosstab")?;

    for (i, path) in [(1usize, &ri_path), (2usize, &ci_path)] {
        output::write_column_header(w, &cols[i], n_cells)?;
        w.i32_list_header(n_cells)?;
        pour(w, path, n_cells * 4)?;
        pagecache::release_path(path);
        let _ = std::fs::remove_file(path);
        rep.at(
            progress::band(progress::WRITE, i, 2),
            "Writing the index columns",
        );
    }
    Ok(())
}

/// Append `i32`s in the little-endian layout that TSON and the spill files share.
fn write_i32_le<W: std::io::Write>(out: &mut W, v: &[i32]) -> Result<()> {
    #[cfg(target_endian = "little")]
    {
        // SAFETY: i32 has no padding, so its bytes are exactly the encoding wanted here.
        let bytes = unsafe { std::slice::from_raw_parts(v.as_ptr() as *const u8, v.len() * 4) };
        out.write_all(bytes)?;
    }
    #[cfg(not(target_endian = "little"))]
    for x in v {
        out.write_all(&x.to_le_bytes())?;
    }
    Ok(())
}

/// Pour a spill file into the result, refusing to if it does not hold what the header promised.
fn pour<W: std::io::Write>(
    w: &mut TsonWriter<W>,
    path: &std::path::Path,
    expect_bytes: usize,
) -> Result<()> {
    use std::io::Read as _;
    let len = std::fs::metadata(path)?.len() as usize;
    if len != expect_bytes {
        return Err(anyhow!(
            "{} holds {len} bytes where the result declares {expect_bytes}",
            path.display()
        ));
    }
    let mut f = std::io::BufReader::with_capacity(1 << 20, std::fs::File::open(path)?);
    let mut buf = vec![0u8; 1 << 20];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        w.raw(&buf[..n])?;
    }
    Ok(())
}

/// Put the cofactors a run used into the task log, so a user sees them without opening a table.
fn report_cofactors(rep: &Reporter, table: &[output::CofactorRow]) {
    for row in table {
        rep.info(format!(
            "{}: cofactor {:.1}{}",
            row.channel,
            row.cofactor,
            if row.unstable { " (unstable)" } else { "" }
        ));
    }
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
