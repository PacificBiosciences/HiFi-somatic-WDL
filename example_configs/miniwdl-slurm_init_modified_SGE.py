# Copyright (c) 2022 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import logging
import os
import shlex
import sys
from contextlib import ExitStack
from typing import Dict, List

from WDL import Type, Value
from WDL.runtime import config
from WDL.runtime.backend.cli_subprocess import _SubprocessScheduler
from WDL.runtime.backend.singularity import SingularityContainer


class SlurmSingularity(SingularityContainer):
    @classmethod
    def global_init(cls, cfg: config.Loader, logger: logging.Logger) -> None:
        # Set resources to maxsize. The base class (_SubProcessScheduler)
        # looks at the resources of the current host, but since we are
        # dealing with a cluster these limits do not apply.
        cls._resource_limits = {
            "cpu": sys.maxsize,
            "mem_bytes": sys.maxsize,
            "time": sys.maxsize,
        }
        _SubprocessScheduler.global_init(cls._resource_limits)
        # Since we run on the cluster, the images need to be placed in a
        # shared directory. The singularity cache itself cannot be shared
        # across nodes, as it can become corrupted when nodes pull the same
        # image. The solution is to pull image to a shared directory on the
        # submit node. If no image_cache is given, simply place a folder in
        # the current working directory.
        if cfg.get("singularity", "image_cache") == "":
            cfg.override(
                {"singularity": {
                    "image_cache": os.path.join(os.getcwd(),
                                                "miniwdl_singularity_cache")
                }}
            )
        SingularityContainer.global_init(cfg, logger)

    @classmethod
    def detect_resource_limits(cls, cfg: config.Loader,
                               logger: logging.Logger) -> Dict[str, int]:
        return cls._resource_limits  # type: ignore

    @property
    def cli_name(self) -> str:
        return "slurm_singularity"

    def process_runtime(self,
                        logger: logging.Logger,
                        runtime_eval: Dict[str, Value.Base]) -> None:
        """Any non-default runtime variables can be parsed here"""
        super().process_runtime(logger, runtime_eval)
        if "time_minutes" in runtime_eval:
            time_minutes = runtime_eval["time_minutes"].coerce(Type.Int()).value
            self.runtime_values["time_minutes"] = time_minutes

        if "slurm_partition" in runtime_eval:
            slurm_partition = runtime_eval["slurm_partition"].coerce(
                Type.String()).value
            self.runtime_values["slurm_partition"] = slurm_partition

        if "gpuCount" in runtime_eval:
            gpuCount = max(1, runtime_eval["gpuCount"].coerce(Type.Int()).value)
            self.runtime_values["gpuCount"] = gpuCount

        if "slurm_constraint" in runtime_eval:
            slurm_constraint = runtime_eval["slurm_constraint"].coerce(
                Type.String()).value
            self.runtime_values["slurm_constraint"] = slurm_constraint

    def _slurm_invocation(self):
        # We use srun as this makes the submitted job behave like a local job.
        # This also gives informative exit codes back, including 253 for out
        # of memory.
        srun_args = [
            "qrsh",
            "-now",
            "no",
            "-N", self.run_id,
        ]

        partition = self.runtime_values.get("slurm_partition", None)
        if partition is not None:
            srun_args.extend(["--partition", partition])

        cpu = self.runtime_values.get("cpu", None)
        if cpu is not None:
            srun_args.extend(["-pe", "def_slot", str(cpu)])

        memory = self.runtime_values.get("memory_reservation", None)
        if memory is not None:
            # Round to the nearest megabyte.
            srun_args.extend(["-l", f"s_vmem={round(memory / (1024 ** 2))}M"])

        gpu = self.runtime_values.get("gpu", None)
        if gpu:
            gpuCount = self.runtime_values.get("gpuCount", 1)
            srun_args.extend(["--gres", f"gpu:{gpuCount}"])

        time_minutes = self.runtime_values.get("time_minutes", None)
        if time_minutes is not None:
            srun_args.extend(["--time", str(time_minutes)])

        slurm_constraint = self.runtime_values.get("slurm_constraint", None)
        if slurm_constraint is not None:
            srun_args.extend(["--constraint", slurm_constraint])

        if self.cfg.has_section("slurm"):
            extra_args = self.cfg.get("slurm", "extra_args")
            if extra_args is not None:
                srun_args.extend(shlex.split(extra_args))
        print(srun_args)
        return srun_args

    def _run_invocation(self, logger: logging.Logger, cleanup: ExitStack,
                        image: str) -> List[str]:
        singularity_command = super()._run_invocation(logger, cleanup, image)

        slurm_invocation = self._slurm_invocation()
        slurm_invocation.extend(singularity_command)
        logger.info("Slurm invocation: " + ' '.join(
            shlex.quote(part) for part in slurm_invocation))
        return slurm_invocation
