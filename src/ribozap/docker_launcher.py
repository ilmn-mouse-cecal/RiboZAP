import subprocess
import logging
from pathlib import Path
from os import path, environ

logger = logging.getLogger(__name__)

def read_mount_file(mount_file: Path):
    try:
        with mount_file.open() as f:
            return [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        logger.error(f"Mount file not found: {mount_file}")
        raise

def get_required_mounts(out_dir: Path, analysis_name: str, other_paths):
    nextflow_logs = out_dir / analysis_name / "logs"
    mounts = [
        f'-v "{out_dir.resolve()}:{out_dir.resolve()}"',
        f'-v "{out_dir.resolve()}/{analysis_name}:/app/nextflow_files"',
        f'-v "{nextflow_logs.resolve()}:{nextflow_logs.resolve()}"',
    ]
    for path in other_paths:
        mounts.append(f'-v {path}:{path}')

    # rRNA index
    rrna_index = Path(environ["HOME"]) / ".rRNA_database_index"
    mounts.append(f'-v "{rrna_index}:/app/idx"')
    return mounts

def build_docker_command(image, docker_options, mount_lines, container_cmd):
    opts = " ".join(docker_options)
    mounts = " ".join(mount_lines)
    return (
        f'docker run --rm {opts} {mounts} '
        f'{image} bash -c "cd /app && {container_cmd}"'
    )

def run_subprocess_command(docker_cmd):
    logger.info("🐳 Running Docker command:")
    logger.info(docker_cmd)
    try:
        subprocess.run(docker_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Docker command failed: {e}")
        raise

def run_docker(
    image: str, 
    mount_file: Path, 
    out_dir: Path, 
    analysis_name: str, 
    container_cmd: str, 
    other_paths,
    cpus=4, 
    memory=16,
    dry_run=False
):
    docker_options = []
    if cpus:
        docker_options.append(f"--cpus={cpus}")
    if memory:
        docker_options.append(f"--memory={memory}g")

    mount_lines = read_mount_file(mount_file)

    # Add required mounts
    mount_lines += get_required_mounts(out_dir, analysis_name, other_paths)

    # Build docker command
    docker_cmd = build_docker_command(image, docker_options, mount_lines, container_cmd)

    # Run docker
    if dry_run:
        return docker_cmd
    else:
        run_subprocess_command(docker_cmd)
