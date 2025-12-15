import tempfile
from pathlib import Path
import pytest
from unittest.mock import patch
from ribozap.docker_launcher import *

def test_read_mount_file(tmp_path):
    mount_file = tmp_path / "mounts.txt"
    content = [
        '-v "/host/data1":"/container/data1"\n',
        '\n',  # blank line should be ignored
        '-v "/host/data2":"/container/data2"\n'
    ]
    mount_file.write_text("".join(content))
    lines = read_mount_file(mount_file)
    assert lines == [
        '-v "/host/data1":"/container/data1"',
        '-v "/host/data2":"/container/data2"'
    ]

def test_read_mount_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_mount_file(tmp_path / "does_not_exist.txt")

def test_get_required_mounts(tmp_path):
    analysis_name = "my_analysis"
    out_dir = tmp_path
    # Set a fake HOME for rRNA index
    mounts = get_required_mounts(out_dir, analysis_name, [])
    assert any(f'-v "{str(out_dir.resolve())}:{str(out_dir.resolve())}"' in m for m in mounts)
    assert any(f'/app/nextflow_files' in m for m in mounts)
    assert any(f'{str(out_dir.resolve())}/{str(analysis_name)}/logs:{str(out_dir.resolve())}/{str(analysis_name)}/logs' in m for m in mounts)
    assert any('/app/idx' in m for m in mounts)

def test_build_docker_command():
    image = "ribozap:latest"
    docker_options = ["--cpus=4", "--memory=8g"]
    mount_lines = ['-v "/a:/b"', '-v "/c:/d"']
    container_cmd = "echo test-command"
    cmd = build_docker_command(image, docker_options, mount_lines, container_cmd)
    assert "docker run --rm" in cmd
    assert "--cpus=4" in cmd
    assert "--memory=8g" in cmd
    assert '-v "/a:/b"' in cmd
    assert '-v "/c:/d"' in cmd
    assert image in cmd
    assert 'bash -c "cd /app && echo test-command"' in cmd

def test_run_subprocess_command_success():
    cmd = "echo test-command"
    with patch("subprocess.run") as mock_run:
        run_subprocess_command(cmd)
        mock_run.assert_called_once_with(cmd, shell=True, check=True)

def test_run_subprocess_command_failure():
    cmd = "exit 1"
    with patch("subprocess.run", side_effect=Exception("fail")) as mock_run:
        with pytest.raises(Exception):
            run_subprocess_command(cmd)
