def get_wc_constraint(key: str = None, add_tmc: bool = False) -> str:
    """Determine wildcard constraints for wildcards that are dependent on
    entries in the config file. For example, the filter wildcard should be
    constrained to the entries under 'filter' in the config file.

    Args:
        key: A key for the config dictionary
        add_tmc: Add 'tmc' as a constraint as well. Useful for too-many-cells portions
            of the pipeline.

    Returns:
        A regex string that defines the wild constraint associated with
        config[key]
    """
    
    if key is None:
        # If key does not exist in the config file, return the regex below instead.
        # This is the default regex: anything that doesn't contain a
        # '\', '.',  '/', or a whitespace
        return "[^\\/\.\s]*"

    elif exists(key, config):
        constraint = "|".join(config[key])
        if add_tmc:
            constraint += "|tmc"

        return f"({constraint})"

    # If key does not exist in the config file, return the regex below instead.
    # This is the default regex: anything that doesn't contain a
    # '\', '.',  '/', or a whitespace
    else:
        return "[^\\/\.\s]*"


def exists(key: str, dictionary: dict) -> bool:
    """Helper function to determine if a key exists in a
    dictionary and isn't empty.

    Args:
        key: A key that is potentially in the dictionary
        dictionary: the dictionary that will be searched

    Returns:
        True if the key exists in the dictionary and the entry is not empty.
    """
    return key in dictionary and dictionary[key]


def extract_wildcards(template: str, target: str) -> dict:
    """Helper function for get_params_instance.

    Creates a dictionary of wildcards where the key is the wildcard name and the value
    is the wildcard's value. The template contains wildcard names and the target
    contains wildcard values. For example, if the template and targets are:

    template = "k~{wildcard1}_alpha~{wildcard2}"
    target = "k~1.5_alpha~2"

    then the returned dictionary is:

    wildcards = {"wildcard1": 1.5, "wildcard2": 2}

    Args:
        template: template string containing wildcard names enclosed between "{}"
        target: target string that follows the template's structure.
    
    Returns:
        Dictionary containing wildcard and wildcard values as key-value pairs.
    """
    # Find all wildcard names in the template.
    names = re.findall("{[^{]*}", template)
    names = [name.strip("{}") for name in names]

    # Find the wildcard values
    q = re.sub("{[^{]*}", "(.*)", template)
    matches = re.search(q, target)

    # Create the output dictionary
    wildcards = dict()
    for i, name in enumerate(names):
        wildcards[name] = matches.group(i + 1)
    
    return wildcards


def get_params_instance(method_name: str, params_str: str, method_dict: dict) -> dict:
    """Get a parameter instance. This is a slightly modified version of Snakemake's 
    implementation of snakemake.utils.Paramspace.instance(). 
    
    Reason for using the modified version:

    Snakemake's implementation could not be used because instance() does not work with 
    the combinatorial features of this pipeline. More specifically, this pipeline uses 
    a dictionary of Paramspaces to store the parameters of each method. 
    The key for each entry is the method_name. method_name is determined from a 
    wildcard. As such, snakemake.utils.Paramspace.instance cannot be used because it 
    doesn't know how to access the dictionary using method_name as a key.

    Args:
        method_name: The name of the method. For example for filtering, a method_name
            can be "IQR". For clustering, it can be "leiden"
        params_str: The parameters for the method as a string. This string is generated 
            using snakemake.utils.Paramspace.instance_patterns
        method_dict: A dictionary where a key is a method_name and the value is the 
            Paramspace associated with that method.
    """

    # If the prameters are noparams or untrackedparams, return an empty dictionary since
    # there are no parameters
    if params_str == "noparams" or params_str == "untrackedparams":
        return {}

    # Get the method's paramspace
    paramspace = method_dict[method_name].paramspace

    # Get the specific parameter instance
    params_dict = extract_wildcards(paramspace.wildcard_pattern, params_str)

    def convert_value_dtype(name, value):
        """Same function as found in snakemake.utils.Paramspace.instance"""
        if paramspace.dataframe.dtypes[name] == bool and value == "False":
            # handle problematic case when boolean False is returned as
            # boolean True because the string "False" is misinterpreted
            return False
        else:
            return pd.Series([value]).astype(paramspace.dataframe.dtypes[name])[0]

    return {
        name: convert_value_dtype(name, value) for name, value in params_dict.items()
    }