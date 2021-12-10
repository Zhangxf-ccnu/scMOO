README file for R package supporting the paper    **"Imputing dropouts for single-cell RNA sequencing based on multi-objective optimization"**.


Contents of this archive
------------------------
This archive contains:   
 
(1) pkg: subdirectory that contains the R package.

(2) scMOO_1.1.pdf: reference manual.

The `scMOO` package has the following R-package dependencies: `keras`, `tensorflow` and `rsvd`.
The dependent packages will be automatically installed along with `scMOO`. You can use the following commands to install `scMOO` from GitHub.

Installation
------------------------
### Step 1. If the devtools package has not been installed, install the devtools package first. Invoke R and then type

`install.packages("devtools")`

### Step 2. Load the devtools package.

`library("devtools")`

### Step 3. Install the scMOO package from GitHub.

`install_github("Zhangxf-ccnu/scMOO", subdir="pkg")`


Quick start
------------------------

### Step 1. Load the library scMOO in R console, by running

`library(scMOO)`

### Step 2. Taking the PBMC_CL dataset as an example, run the following code:

```
data("PBMC_CL")

result <- scMOO(PBMC_CL, percent=0)
```

For detailed usages, please refer to "scMOO_1.1.pdf".

Tutorial
------------------------
A tutorial with examples of cell clustering (downstream experiment) illustrating the usage of `scMOO` is available at:
[scMOO-tutorial.html](https://github.com/Zhangxf-ccnu/scMOO/blob/main/pkg/vignettes/scMOO-tutorial.html) 


Contact
------------------------
Please do not hesitate to contact Miss **Ke Jin** [kej13@mails.ccnu.edu.cn](kej13@mails.ccnu.edu.cn) or Dr. **Xiao-Fei Zhang** [zhangxf@mail.ccnu.edu.cn](zhangxf@mail.ccnu.edu.cn) to seek any clarifications regarding any contents or operation of the archive.
