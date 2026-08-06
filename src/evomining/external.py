"""
evomining.external
==================

Resolution of the external binaries the pipeline shells out to.

The design goal is that in a correctly built `evomining` conda environment,
every tool is simply on $PATH and nothing needs configuring. The old
hard-coded prefixes (/vol/local/conda_envs/muscle, .../fast_tree, ...) are
gone from the code entirely.

If a tool genuinely has to live in a separate environment (an old BLAST for
reproducing a previous run, say), record it in the manifest:

    evomining config --set-tool blastp=/vol/local/conda_envs/blast+

which makes the command run as `micromamba run -p <prefix> blastp ...`.
"""

from __future__ import annotations

import shutil
import subprocess
import sys

# tool -> conda package that provides it, used only for error messages
PROVIDED_BY = {
    "makeblastdb": "blast",
    "blastp":      "blast",
    "muscle":      "muscle",
    "mafft":       "mafft",
    "trimal":      "trimal",
    "FastTree":    "fasttree",
    "fasttree":    "fasttree",
    "nw_clade":    "newick_utils",
    "nw_labels":   "newick_utils",
    "nw_reroot":   "newick_utils",
    "nw_rename":   "newick_utils",
    "nw_display":  "newick_utils",
}


class ToolError(RuntimeError):
    pass


class Tools:
    """Builds argv lists for external tools, honouring per-tool env overrides."""

    def __init__(self, tool_env=None):
        self.tool_env = dict(tool_env or {})
        self._checked = {}

    def argv(self, tool, args):
        """Return the full argv for `tool` with `args`."""
        prefix = self.tool_env.get(tool)
        if prefix:
            return ["micromamba", "run", "-p", str(prefix), tool, *map(str, args)]
        return [tool, *map(str, args)]

    def available(self, tool):
        """True if the tool can be run (on $PATH, or has an env override)."""
        if tool in self._checked:
            return self._checked[tool]
        ok = bool(self.tool_env.get(tool)) or shutil.which(tool) is not None
        self._checked[tool] = ok
        return ok

    def require(self, *tools):
        """Fail early, with an actionable message, if any tool is missing."""
        missing = [t for t in tools if not self.available(t)]
        if missing:
            pkgs = sorted({PROVIDED_BY.get(t, t) for t in missing})
            raise ToolError(
                "required tool(s) not found on $PATH: " + ", ".join(missing) + "\n"
                "  Install them into the active environment:\n"
                "    micromamba install -c conda-forge -c bioconda " + " ".join(pkgs) + "\n"
                "  or point EvoMining at an existing environment, e.g.\n"
                "    evomining config --set-tool " + missing[0] + "=/path/to/conda/env"
            )

    def run(self, tool, args, capture=True, stdin=None, check=False, quiet=False):
        """Run a tool. Returns the CompletedProcess."""
        argv = self.argv(tool, args)
        result = subprocess.run(
            argv,
            input=stdin,
            capture_output=capture,
            text=True,
        )
        if result.returncode != 0 and not quiet:
            err = (result.stderr or "").strip().splitlines()
            tail = err[-1][:200] if err else f"exit {result.returncode}"
            print(f"  WARNING: {tool} failed: {tail}", file=sys.stderr)
        if check and result.returncode != 0:
            raise ToolError(f"{tool} failed (exit {result.returncode}):\n"
                            f"  {' '.join(argv)}\n{result.stderr}")
        return result
