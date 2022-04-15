# List of functions that parse the config files in config/

import copy
import itertools
import os
import pandas as pd
import pathlib

from pathlib import Path
from snakemake.utils import Paramspace
from typing import Final, List, Union, Tuple


class Sample:
    def __init__(
        self, name: str, data_path: str, label_path: str, obs_tag: str
    ) -> None:

        self.name = name
        self.data_path = Path(data_path)
        self.label_path = Path(label_path)
        self.obs_tag = obs_tag

    @classmethod
    def get_all_samples(cls, config: dict) -> List["__class__"]:

        samples = []

        for key in config["inputs"]:
            sample = config["inputs"][key]
            samples.append(
                cls(key, sample["data_path"], sample["label_path"], sample["obs_tag"])
            )

        return samples


class Method:
    def __init__(self, method_type: str, name: str) -> None:
        """Init method for Method class"""

        self.method_type = method_type
        self.name = Method._get_method(method_type, name)
        self.script = Method._get_script(method_type, name)
        self.paramspace = Method._get_paramspace(method_type, name)

    @staticmethod
    def _get_method(method_type: str, method_name: str) -> str:
        is_implemented = method_name in config["implemented"][method_type]
        return method_name if is_implemented else Path(method_name).stem

    @staticmethod
    def _get_script(method_type: str, method_name: str) -> str:
        is_implemented = method_name in config["implemented"][method_type]
        return f"{method_type}.py" if is_implemented else method_name

    @staticmethod
    def _get_paramspace(method_type: str, method_name: str) -> Paramspace:
        is_implemented = method_name in config["implemented"][method_type]

        if config[method_type][method_name]:
            params = copy.deepcopy(config[method_type][method_name])

        elif is_implemented and config["implemented"][method_type][method_name]:
            params = copy.deepcopy(config["implemented"][method_type][method_name])

        else:
            params = {}

        for key, val in params.items():
            if not isinstance(val, list):
                params[key] = [val]

        # Get all possible parameter combinations and store them in a list
        all_combs = [list(i) for i in zip(*itertools.product(*params.values()))]
        df_all_combs = pd.DataFrame(dict(zip(params.keys(), all_combs)))

        return Paramspace(df_all_combs, filename_params="*")

    @classmethod
    def get_methods(cls, method_type: str) -> List["__class__"]:
        return [cls(method_type, method_name) for method_name in config[method_type]]


class Tree:
    def __init__(self, data: object = None) -> None:
        self.data = data
        self.children = []

    @staticmethod
    def get_paths(
        t: "__class__",
        paths: List["__class__"] = None,
        current_path: List["__class__"] = None,
    ) -> List["__class__"]:

        if paths is None:
            paths = []

        if current_path is None:
            current_path = []

        if t.data is not None:
            current_path.append(t.data)

        if len(t.children) == 0:
            paths.append(current_path)

        else:
            for child in t.children:
                Tree.get_paths(child, paths, list(current_path))

        return paths


def make_method_tree(methodorder: List[str] = []) -> Tree:

    method_tree = Tree()

    def helper(tree: Tree, methodorder):
        if len(methodorder) == 0:
            return tree

        method_type = methodorder[0]

        if not exists(method_type, config):
            return tree

        for method_name in config[method_type]:
            method = Method(method_type, method_name)

            if len([*method.paramspace.instance_patterns]) == 0:
                tree.children.append(
                    Tree(data=f"{method.method_type}_{method.name}_noparams")
                )
            else:
                for param_instance in method.paramspace.instance_patterns:
                    tree.children.append(
                        Tree(
                            data=(
                                f"{method.method_type}_{method.name}"
                                f"_{param_instance}"
                            )
                        )
                    )

        for child in tree.children:
            helper(child, methodorder[1:])

    helper(method_tree, methodorder)

    return method_tree


def get_samples_basenames(skip_merged: bool = False) -> List[str]:

    basenames = []
    samples = [sample.name for sample in Sample.get_all_samples(config)]

    if len(samples) > 1:
        if config["run_on_each_sample"]:
            basenames.extend(
                expand(f"{config['output_dir']}/{{sample}}", sample=samples)
            )

        if not skip_merged:
            basenames.append(f"{config['output_dir']}/merged")

    elif len(samples) == 1:
        basenames.append(f"{config['output_dir']}/{samples[0]}")

    else:
        raise ValueError("Minimum of one input required.")

    return [Path(b) for b in basenames]


def parse_methods(
    methodorder: List[str] = [], skip_merged: bool = False, filename: str = "mtx.h5ad"
) -> List[str]:

    if len(methodorder) == 0:
        return []
    elif not exists(methodorder[0], config):
        return []

    output_files = []
    first_method_type = methodorder[0]
    method_tree = make_method_tree(methodorder=methodorder)

    for parent_dir in get_samples_basenames(skip_merged=skip_merged):
        for method_param_path in Tree.get_paths(method_tree):
            path = "/".join(method_param_path)

            output_files.append(parent_dir / first_method_type / path / filename)

    return output_files


def get_plots():
    """Determine the list of desired plots based off what the user
    defined in the config files

    Returns:
        List of strings where each element is a plot file
    """
    if not exists("plot", config):
        return []

    output_files = []
    use_labels = config["plot"]["use_labels"]

    if config["plot"]["dim_reduce"] and exists("dim_reduce", config):

        method_tree = make_method_tree(methodorder=["dim_reduce", "filter", "norm"])
        for parent_dir in get_samples_basenames():
            for method_param_path in Tree.get_paths(method_tree):
                path = "/".join(method_param_path)

                if use_labels:
                    for label in config["labels"]:
                        output_files.append(
                            parent_dir
                            / "figures"
                            / "dim_reduce"
                            / path
                            / f"scatter_{label}.html"
                        )

                output_files.append(
                    parent_dir / "figures" / "dim_reduce" / path / "scatter.html"
                )

    if config["plot"]["cluster"] and exists("cluster", config):

        method_tree = make_method_tree(
            methodorder=["cluster", "filter", "norm", "dim_reduce"]
        )
        for parent_dir in get_samples_basenames():
            for method_param_path in Tree.get_paths(method_tree):
                path = "/".join(method_param_path)
                output_files.append(
                    parent_dir / "figures" / "cluster" / path / "scatter.html"
                )

    return output_files


def get_too_many_cells_output():
    if not exists("too-many-cells", config):
        return []

    output_files = []
    dict_tmc = config["too-many-cells"]

    for key in dict_tmc:
        if "make-tree" in key:
            # If make-tree uses a prior
            if exists("prior", dict_tmc[key]):
                if dict_tmc[dict_tmc[key]["prior"]]["tmc_qc"]:
                    for parent_dir in get_samples_basenames():
                        output_files.append(
                            parent_dir
                            / "too-many-cells"
                            / key
                            / "filter_tmc_untrackedparams"
                            / "norm_tmc_untrackedparams"
                            / f"{dict_tmc[key]['prior']}.prior.{key}.done"
                        )

                else:
                    method_tree = make_method_tree(methodorder=["filter", "norm"])
                    for parent_dir in get_samples_basenames():
                        for method_param_path in Tree.get_paths(method_tree):
                            path = "/".join(method_param_path)
                            output_files.append(
                                parent_dir
                                / "too-many-cells"
                                / key
                                / path
                                / f"{dict_tmc[key]['prior']}.prior.{key}.done"
                            )

            # If the make-tree doesn't contain a prior
            else:
                if dict_tmc[key]["tmc_qc"]:
                    for parent_dir in get_samples_basenames():
                        output_files.append(
                            parent_dir
                            / "too-many-cells"
                            / key
                            / "filter_tmc_untrackedparams"
                            / "norm_tmc_untrackedparams"
                            / f"{key}.done"
                        )
                else:
                    method_tree = make_method_tree(methodorder=["filter", "norm"])
                    for parent_dir in get_samples_basenames():
                        for method_param_path in Tree.get_paths(method_tree):
                            path = "/".join(method_param_path)
                            output_files.append(
                                parent_dir
                                / "too-many-cells"
                                / key
                                / path
                                / f"{key}.done"
                            )

    return output_files


def get_differential_output() -> List[str]:
    if not exists("differential", config):
        return []

    output_files = []

    method_tree = make_method_tree(
        methodorder=["filter", "norm", "dim_reduce", "cluster"]
    )
    for parent_dir in get_samples_basenames():
        for method_param_path in Tree.get_paths(method_tree):
            path = "/".join(method_param_path)
            for diff in [*config["differential"].keys()]:
                output_files.append(
                    parent_dir
                    / "differential"
                    / f"differential_{diff}"
                    / path
                    / "mtx.h5ad"
                )

    return output_files


def get_final_output() -> List[str]:

    output_files = []

    methodorder = ["filter", "norm", "dim_reduce", "cluster"]
    for i, method_type in reversed(list(enumerate(methodorder))):

        methodorder_shuffled = [method_type] + methodorder[:i]

        if method_type == "filter" and exists(method_type, config):
            output_files.extend(
                parse_methods(methodorder=methodorder_shuffled, skip_merged=True)
            )

        elif exists(method_type, config):
            output_files.extend(parse_methods(methodorder=methodorder_shuffled))
            break

    output_files.extend(get_plots())
    output_files.extend(get_too_many_cells_output())
    output_files.extend(get_differential_output())

    return output_files


def get_method_dict(method_type: str) -> dict:
    if not exists(method_type, config):
        return {}

    method_dict = dict()

    methods = Method.get_methods(method_type)
    for method in methods:
        method_dict[method.name] = method

    return method_dict


def get_differential_script(wildcards) -> str:
    diff = config['differential'][wildcards.diff]

    if exists("script", diff):
        return diff['script']
    else:
        return "differential.py"
    
