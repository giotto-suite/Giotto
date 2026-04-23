# How to Contribute?

We welcome contributions or suggestions from other developers. Please contact us if you have questions or would like to discuss an addition or major modifications to the Giotto main code.
The source code for Giotto Suite may be found on our [GitHub repository](https://github.com/giotto-suite/Giotto/).

$~$

### Coding Style  
Following a particular programming style will help programmers read and understand source code conforming to the style, and help to avoid introducing errors. Here we present a small list of guidelines on what is considered a good practice when writing R codes in Giotto package. Most of them are adapted from [Bioconductor - coding style](https://bioconductor.org/developers/how-to/coding-style/) or [Google's R Style Guide](https://google.github.io/styleguide/Rguide.xml). These guidelines are preferences and strongly encouraged!  

- **Naming**
    - Use camelCase for Giotto user-facing exported function names. (`functionName()`)  
    - Use snake_case for non-user-facing exported functions, which are essentially any functions not directly related to commonly used data processing, analysis, and visualization. (`function_name()`)  
    - Use "." prefix and snake_case for internal non-exported functions. (`.function_name()`)  
    - Use snake_case for parameter names.  
    - Do not use "." as a separator in function naming. (in the S3 class system, some(x) where x is class A will dispatch to some.A)  

- **Use of space**
    - Do not place a space before a comma, but always place one after a comma. This: `a, b, c`.  
    - Always use space around “=” when using named arguments to functions. This: `somefunc(a = 1, b = 2)`.
    
- **Use of symbols**
    - Do not use any non-UTF-8 characters unless provided as the escape code. For example: `'\u00F6'` for `ö`
    
Beyond these guidelines, [*styler*](https://github.com/r-lib/styler) should be used in order to maintain code uniformity.

$~$

### Stat functions  
Most Giotto commands can accept several matrix classes (DelayedMatrix, SparseM, Matrix or base matrix). To 
facilitate this we provide **flexible** wrappers that work on any type of matrix class.

- **mean_flex**:        analogous to mean()  
- **rowSums_flex**:     analogous to rowSums()  
- **rowMeans_flex**:    analogous to rowMeans()  
- **colSums_flex**:     analogous to colSums()  
- **colMeans_flex**:    analogous to colMeans()  
- **t_flex**:           analogous to t()  
- **cor_flex**:         analogous to cor()  

$~$

### Auxiliary functions  
Giotto has a number of auxiliary or convenience functions that might help you to adapt your code or write new code for Giotto. We encourage you to use these small functions to maintain uniformity throughout the code.

- **lapply_flex**: analogous to lapply() and works for both windows and unix systems  
- **all_plots_save_function**: compatible with Giotto instructions and helps to automatically save generated plots  
- **plot_output_handler**: further wraps `all_plots_save_function` and includes handling for return_plot and show_plot and Giotto instructions checking  
- **determine_cores**: to determine the number of cores to use if a user does not set this explicitly  
- **get_os**: to identify the operating system  
- **update_giotto_params**: will catch and store the parameters for each used command on a giotto object  
- **wrap_txt** and **wrap_msg**: text and message formatting functions  
- **vmsg**: framework for Giotto's verbosity-flagged messages  
- **package_check**: to check if a package exists, works for packages on CRAN, Bioconductor and Github

The last function should be used within your contribution code. It has the additional benefit that it will suggest the user how to download the package if it is not available. To keep the size of Giotto within limits we prefer not to add too many new dependencies.

$~$

### Package Imports
Giotto tracks packages and functions to import in a centralized manner. When adding code that requires functions from another package, add the roxygen tags to the `package_imports.R` file for that Giotto module.

$~$

### Getters and Setters
Giotto stores information in different [slots](articles/structure.html#giotto-object-structure), which can be accessed through these getters and setters functions. They can be found in the [accessors.R](https://github.com/giotto-suite/Giotto/blob/suite/R/accessors.R) file. 

- **getCellMetadata()**: Gets cell metadata  
- **setCellMetadata()**: Sets cell metadata  

- **getFeatureMetadata()**: Gets feature metadata  
- **getFeatureMetadata()**: Sets feature metadata  

- **getExpression()**: To select the expression matrix to use  
- **setExpression()**: Sets a new expression matrix to the expression slot  

- **getSpatialLocations()**: Get spatial locations to use  
- **setSpatialLocations()**: Sets new spatial locations  

- **getDimReduction()**: To select the dimension reduction values to use  
- **setDimReduction()**: Sets new dimension reduction object  

- **getNearestNetwork()**: To select the nearest neighbor network (kNN or sNN) to use  
- **setNearestNetwork()**: Sets a new nearest neighbor network (kNN or sNN)  

- **getSpatialNetwork()**: To select the spatial network to use  
- **setSpatialNetwork()**: Sets a new spatial network  

- **getPolygonInfo()**: Gets spatial polygon information  
- **setPolygonInfo()**: Set new spatial polygon information  

- **getFeatureInfo()**: Gets spatial feature information  
- **setFeatureInfo()**: Sets new spatial feature information  

- **getSpatialEnrichment()**: Gets spatial enrichment information  
- **setSpatialEnrichment()**: Sets new spatial enrichment information  

- **getMultiomics()**: Gets multiomics information  
- **setMultiomics()**: Sets multiomics information  

$~$

### Python code
To use Python code we prefer to create a python wrapper/functions around the python code, which can then be sourced by reticulate. As an example we show the basic principles of how we implemented the Leiden clustering algorithm.  

1. write python wrapper and store as python_leiden.py in */inst/python*:    
```
import igraph as ig
import leidenalg as la
import pandas as pd
import networkx as nx

def python_leiden(df, partition_type, initial_membership=None, weights=None, n_iterations=2, seed=None, resolution_parameter = 1):
    
    # create networkx object
    Gx = nx.from_pandas_edgelist(df = df, source = 'from', target =  'to', edge_attr = 'weight')  
    
    # get weight attribute
    myweights = nx.get_edge_attributes(Gx, 'weight')

    ....

    return(leiden_dfr)
```


2. source python code with reticulate:  
```
python_leiden_function = system.file("python", "python_leiden.py", package = 'Giotto')
reticulate::source_python(file = python_leiden_function)
```

3. use python code as if R code:  
See **doLeidenCLuster** for more detailed information.  
```
 pyth_leid_result = python_leiden(df = network_edge_dt,
                                   partition_type = partition_type,
                                   initial_membership = init_membership,
                                   weights = 'weight',
                                   n_iterations = n_iterations,
                                   seed = seed_number,
                                   resolution_parameter = resolution)
```

