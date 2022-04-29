import copy
import itertools
import os
import pandas as pd
import pathlib

from pathlib import Path
from snakemake.utils import Paramspace
from typing import Final, List, Union, Tuple


class Sample:
    """Holds all the necessary information for a sample. A sample is a single entry
    under config["inputs"] (e.g input_1)

    Attributes:
        name: The name of the sample (e.g input_1, input_2, etc.)
        data_path: Path to directory containing the feature barcode matrix for that
            sample
        label_path: Path to the label file. A label file is a csv file that contains
            annotations for the observation. Only categorical/binary variables are
            allowed.
        obs_tag: A tag that will be appended to the end of each barcode.
    """
    def __init__(
        self, name: str, data_path: str, label_path: str, obs_tag: str
    ) -> None:
        self.name = name
        self.data_path = Path(data_path)
        self.label_path = Path(label_path)
        self.obs_tag = obs_tag

    @classmethod
    def get_all_samples(cls, config: dict) -> List["__class__"]:
        """Get all the samples provided in the config file.
        
        Args:
            config: config file in dictionary format.
        
        Returns:
            A list of Sample objects
        """

        samples = []

        for key in config["inputs"]:
            sample = config["inputs"][key]
            samples.append(
                cls(key, sample["data_path"], sample["label_path"], sample["obs_tag"])
            )

        return samples


class Method:
    """Stores all the important information about a method.
    
    Attributes:
        method_type: The type of method (e.g, filter, norm, cluster, etc.)
        name: The name of the method (e.g IQR, leiden, etc.)
        script: The name of the script that will run the method
        paramspace: An instance of Paramspace for that particular method.
    """

    def __init__(self, method_type: str, name: str) -> None:
        self.method_type = method_type
        self.name = Method._get_method(method_type, name)
        self.script = Method._get_script(method_type, name)
        self.paramspace = Method._get_paramspace(method_type, name)

    @staticmethod
    def _get_method(method_type: str, method_name: str) -> str:
        """Helper method that determines the method from the config file.
        
        Args:
            method_type: The type of method (e.g, filter, norm, cluster, etc.)
            method_name: The name of the method (e.g IQR, leiden, etc.)
        
        Returns:
            A string representing the method name.
        """
        is_implemented = method_name in config["implemented"][method_type]

        # If the method_name is a stript, return only the stem
        return method_name if is_implemented else Path(method_name).stem

    @staticmethod
    def _get_script(method_type: str, method_name: str) -> str:
        """Helper method that determines the method script from the config file.
        
        Args:
            method_type: The type of method (e.g, filter, norm, cluster, etc.)
            method_name: The name of the method (e.g IQR, leiden, etc.)
        
        Returns:
            A string representing the name of the script.
        """
        is_implemented = method_name in config["implemented"][method_type]
        return f"{method_type}.py" if is_implemented else method_name

    @staticmethod
    """Helper method that determines the method Paramspace from the config file.
        
        Args:
            method_type: The type of method (e.g, filter, norm, cluster, etc.)
            method_name: The name of the method (e.g IQR, leiden, etc.)
        
        Returns:
            An instance of Paramspace.
        """
    def _get_paramspace(method_type: str, method_name: str) -> Paramspace:
        is_implemented = method_name in config["implemented"][method_type]

        # If the method exists create a copy of the dictionary of parameters
        # from the config file. Otherwise, set params to be an empty dictionary
        if config[method_type][method_name]:
            params = copy.deepcopy(config[method_type][method_name])

        elif is_implemented and config["implemented"][method_type][method_name]:
            params = copy.deepcopy(config["implemented"][method_type][method_name])

        else:
            params = {}

        # For any parameters that are not applied combinatorially, wrap the parameters
        # in a one-element list. This is to ensure that all values in the params
        # dictionary is a list.
        for key, val in params.items():
            if not isinstance(val, list):
                params[key] = [val]

        # Get all possible parameter combinations and store them in a list
        all_combs = [list(i) for i in zip(*itertools.product(*params.values()))]
        df_all_combs = pd.DataFrame(dict(zip(params.keys(), all_combs)))

        return Paramspace(df_all_combs, filename_params="*")

    @classmethod
    def get_methods(cls, method_type: str) -> List["__class__"]:
        """Retrieve all the methods for a particular method_type.
        
        Args:
            method_type: The type of method (e.g, filter, norm, cluster, etc.)

        Returns:
            A list of Methods
        """
        return [cls(method_type, method_name) for method_name in config[method_type]]


class Tree:
    """Helper class that implements a basic tree data structure. This tree
    implementation is defined recursively. In other words, a node is also a Tree.

    Attributes:
        data: The data stored at the node.
        children: A list of children nodes. These nodes are instances of Tree.
    """

    def __init__(self, data: object = None) -> None:
        self.data = data
        self.children = []

    @staticmethod
    def get_paths(
        t: "__class__",
        paths: List["__class__"] = None,
        current_path: List["__class__"] = None,
    ) -> List["__class__"]:
        """Get all possible paths of the tree starting from the root node and ending
        at the leaves. This implementation is recursive.

        Args:
            t: An instance of Tree.
            paths: An accumulating list of completed paths.
            current_path: The current path being searched.
        """

        if paths is None:
            paths = []

        if current_path is None:
            current_path = []

        # If the data is not None, then append the data to the current path
        if t.data is not None:
            current_path.append(t.data)

        # If the current node has no children, then we have reached a leaf node and 
        # thus the current path is completed. Append this path to paths.
        if len(t.children) == 0:
            paths.append(current_path)

        # If the current node has children, then get the paths from the current node
        # to the children.
        else:
            for child in t.children:
                Tree.get_paths(child, paths, list(current_path))

        return paths


def make_method_tree(methodorder: List[str] = []) -> Tree:
    """Create a tree of methods and their parameters. This function builds the output
    directory structure. Each node is the name of the folder. The name of the folder
    follows the following structure:
        {method_type}_{method_name}_{parameter_instance}
    
    e.g
        filter_IQR_k~1.5
    
    Args:
        methodorder: A list of method_types that will the determine which method_types
            will be found higher up on the list. Method_types at the beginnig of the 
            list will be found closer to the root node. The root node is the first
            method_type in the list.
    
    Returns:
        A Tree of Methods. Each node stores a Method

    """

    method_tree = Tree()

    def helper(tree: Tree, methodorder) -> None:
        """Recursive helper function that actually creates the method tree.
        make_method_tree is just a wrapper.

        Args:
            tree: The tree being built.
            methodorder: The order in which the methods will be found in the tree.

        """

        # If there is no methoorder, then don't do anything
        if len(methodorder) == 0:
            return

        method_type = methodorder[0]

        # If the method_type does not exist, then don't proceed any further. The
        # method tree that was built so far will be returned (by the outer function)
        if not exists(method_type, config):
            return

        # Go through each method type and append each possible combination of method
        # and their parameters to the tree.
        for method_name in config[method_type]:
            method = Method(method_type, method_name)

            # If there are no parameters, then indicate "no_params"
            if len([*method.paramspace.instance_patterns]) == 0:
                tree.children.append(
                    Tree(data=f"{method.method_type}_{method.name}_noparams")
                )

            # Otherwise, determine all possible parameter combination and append it to 
            # the tree
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

        # Run this function on the next method_type in methodorder
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


def get_method_dict(method_type: str) -> dict:
    """Get a dictionary containing all the methods for a specific method_type.
    A method_type for example can be "filter" and a method for filter would be "IQR"
    """

    # If the method_type does not exist, return an empty dictionary
    if not exists(method_type, config):
        return {}

    # Loop through the list of methods and create a dictionary. Each key is a method
    # name, and each value is an Method instance.
    method_dict = dict()
    methods = Method.get_methods(method_type)
    for method in methods:
        method_dict[method.name] = method

    return method_dict


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
                if exists("tmc_qc", dict_tmc[dict_tmc[key]["prior"]]):
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
                diff_dict = config["differential"][diff]
                
                if "one_vs_rest" in diff_dict:
                    if diff_dict["one_vs_rest"] in [*config["cluster"].keys()]
                        if diff_dict["one_vs_rest"] in path:                
                            output_files.append(
                                parent_dir
                                / "differential"
                                / f"differential_{diff}"
                                / path
                                / "mtx.h5ad"
                            )
                        else:
                            continue
                    else:
                        if "filters" in diff_dict["one_vs_rest"]:
                            s1 = set([*diff_dict["one_vs_rest"]["filters"].keys()])
                            set_cluster = set([*config["cluster"].keys])
                            if len(s1.intersection(set_cluster)) != 0:
                                cluster_method = list(s1.intersection(set_cluster))[0]
                                if cluster_method is path:
                                    output_files.append(
                                        parent_dir
                                        / "differential"
                                        / f"differential_{diff}"
                                        / path
                                        / "mtx.h5ad"
                                    )
                                else:
                                    continue 
                elif "group1" in diff_dict and "group2" in diff_dict:
                    s1 = set([*diff_dict["group1"].keys()])
                    s2 = set([*diff_dict["group2"].keys()])
                    set_cluster = set([*config["cluster"].keys])                   
                    
                    if len(s1.intersection(s2, set_cluster)) != 0:
                        cluster_method = list(s1.intersection(s2, set_cluster))[0]
                        if cluster_method is path:
                            output_files.append(
                                parent_dir
                                / "differential"
                                / f"differential_{diff}"
                                / path
                                / "mtx.h5ad"
                            )
                        else:
                            continue

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
    
