//! Where does the time go when reading a crosstab? Fetch and decode, timed separately.
//!
//!     TERCEN_URI=… TERCEN_TOKEN=… probe_read <tableId> <cells> [chunk]
use std::sync::Arc;
use std::time::Instant;

use anyhow::Result;
use tercen_rs::TercenClient;

#[tokio::main]
async fn main() -> Result<()> {
    asinh_operator::init_tracing();
    let a: Vec<String> = std::env::args().collect();
    let (table, cells) = (a[1].clone(), a[2].parse::<usize>()?);
    let chunk: usize = a.get(3).and_then(|v| v.parse().ok()).unwrap_or(1_000_000);
    let client = Arc::new(
        TercenClient::from_env()
            .await
            .map_err(|e| anyhow::anyhow!("{e}"))?,
    );
    let streamer = tercen_rs::table::TableStreamer::new(&client);

    let (mut fetch, mut decode, mut got, mut bytes) = (0.0, 0.0, 0usize, 0usize);
    let t0 = Instant::now();
    let mut offset = 0usize;
    while offset < cells {
        let want = chunk.min(cells - offset);
        let t = Instant::now();
        let raw = streamer
            .stream_tson(
                &table,
                Some(vec![".ri".into(), ".ci".into(), ".y".into()]),
                offset as i64,
                want as i64,
            )
            .await
            .map_err(|e| anyhow::anyhow!("{e}"))?;
        fetch += t.elapsed().as_secs_f64();
        bytes += raw.len();
        let t = Instant::now();
        let c = asinh_operator::input::decode_chunk(&raw, true, true, true)?;
        decode += t.elapsed().as_secs_f64();
        println!(
            "  chunk: asked {want}, got {} rows, {:.1} MB ({:.0} B/row)",
            c.len(),
            raw.len() as f64 / 1e6,
            raw.len() as f64 / c.len().max(1) as f64
        );
        if c.is_empty() {
            break;
        }
        got += c.len();
        offset += c.len();
    }
    let total = t0.elapsed().as_secs_f64();
    println!(
        "cells {got} in {total:.1}s  fetch {fetch:.1}s  decode {decode:.1}s  \
         ({:.0} cells/s, {:.1} MB over the wire)",
        got as f64 / total,
        bytes as f64 / 1e6
    );
    Ok(())
}
