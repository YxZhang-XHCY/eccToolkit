"""Subprocess utilities for eccToolkit."""

import logging
import shutil
import subprocess
from typing import List, Optional, Union

logger = logging.getLogger(__name__)


def run_command(
    cmd: Union[str, List[str]],
    shell: bool = True,
    capture_output: bool = True,
    check: bool = False,
    timeout: Optional[int] = None,
    cwd: Optional[str] = None,
) -> subprocess.CompletedProcess:
    """
    Run an external command.

    Args:
        cmd: Command string or list
        shell: Run through shell
        capture_output: Capture stdout/stderr
        check: Raise CalledProcessError on failure
        timeout: Timeout in seconds
        cwd: Working directory

    Returns:
        CompletedProcess instance
    """
    logger.debug(f"Running command: {cmd}")

    try:
        result = subprocess.run(
            cmd,
            shell=shell,
            capture_output=capture_output,
            text=True,
            timeout=timeout,
            cwd=cwd,
        )

        if result.returncode != 0:
            logger.warning(f"Command returned {result.returncode}")
            if result.stderr:
                logger.warning(f"stderr: {result.stderr.strip()}")

        if check and result.returncode != 0:
            raise subprocess.CalledProcessError(
                result.returncode, cmd, result.stdout, result.stderr
            )

        return result

    except subprocess.TimeoutExpired:
        logger.error(f"Command timed out after {timeout}s: {cmd}")
        raise
    except Exception as e:
        logger.error(f"Command failed: {e}")
        raise


def check_tool_installed(tool_name: str) -> bool:
    """
    Check if a command-line tool is installed.

    Args:
        tool_name: Name of the tool

    Returns:
        True if tool is available
    """
    return shutil.which(tool_name) is not None


def get_tool_version(tool_name: str, version_flag: str = "--version") -> Optional[str]:
    """
    Get version string for a tool.

    Args:
        tool_name: Name of the tool
        version_flag: Flag to get version

    Returns:
        Version string or None
    """
    if not check_tool_installed(tool_name):
        return None

    try:
        result = run_command(f"{tool_name} {version_flag}", check=False)
        if result.returncode == 0:
            return result.stdout.strip().split("\n")[0]
    except Exception:
        pass

    return None


def require_tools(tools: List[str]) -> None:
    """
    Check that required tools are installed.

    Args:
        tools: List of tool names

    Raises:
        RuntimeError: If any tool is missing
    """
    missing = []
    for tool in tools:
        if not check_tool_installed(tool):
            missing.append(tool)

    if missing:
        raise RuntimeError(
            f"Missing required tools: {', '.join(missing)}. "
            "Please install them and ensure they are in your PATH."
        )

    logger.info(f"All required tools available: {', '.join(tools)}")
