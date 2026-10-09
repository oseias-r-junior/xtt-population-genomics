#!/usr/bin/env python3
import argparse
import subprocess
import sys
import time
import shutil
import socket
from pathlib import Path
from datetime import datetime

EXPECTED_OUTPUTS: tuple[str, ...] = (
    "contig_report.txt",
    "biomarkers.blast.txt",
    "mge.report.txt",
    "chromosome.fasta",
)

def iso():
    return datetime.now().isoformat(timespec="seconds")

def log(msg):
    print(msg, flush=True)

def check_mob_recon():
    """Verifica se mob_recon está disponível e imprime versão + caminho."""
    path = shutil.which("mob_recon")
    if path is None:
        log("ERROR: mob_recon não encontrado no PATH.")
        sys.exit(1)

    log(f"mob_recon binary: {path}")

    try:
        version = subprocess.check_output(
            ["mob_recon", "--version"],
            text=True,
            encoding="utf-8",
            errors="replace"
        ).strip()
        log(f"mob_recon version: {version}")
    except Exception as e:
        log(f"ERROR ao obter versão do mob_recon: {e}")
        sys.exit(1)

def is_completed(outdir: Path):
    """Verifica se todos os arquivos esperados existem e têm tamanho > 0."""
    for fname in EXPECTED_OUTPUTS:
        f = outdir / fname
        if not f.exists() or f.stat().st_size == 0:
            return False
    return True

def run_command(cmd, logfile, t0):
    """Executa o comando e registra stdout/stderr no log."""
    with open(logfile, "w") as logf:
        logf.write(f"COMMAND: {' '.join(cmd)}\n")
        logf.write(f"HOST: {socket.gethostname()}\n")
        logf.write(f"CWD: {Path.cwd()}\n")
        logf.write(f"Python: {sys.version.split()[0]}\n")
        logf.write(f"Started: {iso()}\n\n")

        try:
            proc = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding="utf-8",
                errors="replace",
                check=False,
            )
        except KeyboardInterrupt:
            raise
        except Exception as e:
            logf.write(f"\nERROR: {e}\n")
            return False, -1

        logf.write("=== STDOUT ===\n")
        logf.write(proc.stdout + "\n")
        logf.write("=== STDERR ===\n")
        logf.write(proc.stderr + "\n")

        elapsed = time.time() - t0
        logf.write(f"\nElapsed: {elapsed:.2f} seconds\n")
        logf.write(f"Finished: {iso()}\n")

        return proc.returncode == 0, proc.returncode

def main():
    parser = argparse.ArgumentParser(description="Batch runner for MOB-SUITE (mob_recon).")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--start-after", type=str, default=None)
    args = parser.parse_args()

    input_dir = Path(args.input)
    output_dir = Path(args.output)
    logs_dir = output_dir.parent / "logs"
    results_dir = output_dir.parent / "results"

    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    if not input_dir.exists():
        log(f"ERROR: input directory not found: {input_dir}")
        sys.exit(1)

    log("========================================")
    log("MOB-SUITE RUN")
    log("========================================")
    log(f"Host: {socket.gethostname()}")
    log(f"Started: {iso()}")
    log(f"Python version: {sys.version.split()[0]}")
    log(f"Threads: {args.threads}")
    log(f"CWD: {Path.cwd()}")
    log("========================================")

    check_mob_recon()

    genomes = sorted(input_dir.glob("*.fna"))
    if len(genomes) == 0:
        log("ERROR: no .fna genomes found.")
        sys.exit(1)

    total = len(genomes)
    log(f"Genomes detected: {total}")

    completed = 0
    skipped = 0
    failed = 0
    failed_list = []
    skipped_list = []

    start_time = time.time()
    processed = 0

    # Normalizar start-after
    target = None
    if args.start_after:
        target = Path(args.start_after).stem.lower()

    started_processing = target is None

    summary_tsv = results_dir / "run_summary.tsv"
    with open(summary_tsv, "w") as tsv:
        tsv.write("genome\tstatus\ttime_seconds\texit_code\n")

    for idx, genome in enumerate(genomes, start=1):
        genome_name = genome.stem
        genome_norm = genome_name.lower()

        # START-AFTER corrigido
        if not started_processing:
            if genome_norm == target:
                started_processing = True
                log(f"[{idx}/{total}] SKIP (start-after match) {genome_name}")
                continue
            else:
                log(f"[{idx}/{total}] SKIP UNTIL MATCH {genome_name}")
                continue

        # LIMIT
        if args.limit is not None and processed >= args.limit:
            log("Limit reached.")
            break

        outdir = output_dir / genome_name
        outdir.mkdir(exist_ok=True)
        logfile = logs_dir / f"{genome_name}.log"

        # SKIP automático
        if is_completed(outdir):
            log(f"[{idx}/{total}] SKIPPED {genome_name} (already completed)")
            skipped += 1
            skipped_list.append(genome_name)
            with open(summary_tsv, "a") as tsv:
                tsv.write(f"{genome_name}\tskipped\t0\t0\n")
            continue

        log(f"[{idx}/{total}] RUN {genome_name}")
        t0 = time.time()

        cmd = [
            "mob_recon",
            "-i", str(genome),
            "-o", str(outdir),
            "--force",
            "-n", str(args.threads),
        ]

        ok, exit_code = run_command(cmd, logfile, t0)
        elapsed = time.time() - t0

        # Checagem de integridade
        if ok and is_completed(outdir):
            log(f"[DONE] {genome_name} ({elapsed:.1f}s)")
            completed += 1
            status = "done"
        else:
            log(f"[FAILED] {genome_name}")
            failed += 1
            failed_list.append(genome_name)
            status = "failed"

        with open(summary_tsv, "a") as tsv:
            tsv.write(f"{genome_name}\t{status}\t{elapsed:.2f}\t{exit_code}\n")

        processed += 1

    # Detectar se start-after nunca apareceu
    if args.start_after and not started_processing:
        log(f"ERROR: genome '{args.start-after}' not found.")
        sys.exit(1)

    total_time = time.time() - start_time

    # ⭐ success_rate corrigido
    success_rate = ((completed + skipped) / total) * 100

    summary_txt = results_dir / "run_summary.txt"
    with open(summary_txt, "w") as sf:
        sf.write("========================================\n")
        sf.write("MOB-SUITE SUMMARY\n")
        sf.write("========================================\n")
        sf.write(f"Genomes found ............. {total}\n")
        sf.write(f"Completed ................. {completed}\n")
        sf.write(f"Skipped ................... {skipped}\n")
        sf.write(f"Failed .................... {failed}\n")
        sf.write(f"Success rate .............. {success_rate:.2f}%\n")
        sf.write(f"Elapsed time .............. {total_time/3600:.2f}h\n")
        sf.write("========================================\n")

        if failed_list:
            sf.write("FAILED GENOMES:\n")
            for g in failed_list:
                sf.write(f" - {g}\n")

        if skipped_list:
            sf.write("\nSKIPPED GENOMES:\n")
            for g in skipped_list:
                sf.write(f" - {g}\n")

    log("\n========================================")
    log("MOB-SUITE SUMMARY")
    log("========================================")
    log(f"Genomes found ............. {total}")
    log(f"Completed ................. {completed}")
    log(f"Skipped ................... {skipped}")
    log(f"Failed .................... {failed}")
    log(f"Success rate .............. {success_rate:.2f}%")
    log(f"Elapsed time .............. {total_time/3600:.2f}h")
    log("========================================")

    if failed > 0:
        sys.exit(1)
    sys.exit(0)

if __name__ == "__main__":
    main()