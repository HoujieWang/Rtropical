# Resubmission 

I got some instructions from CRAN and fixed them as below. 

## 1  
```
If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking.
(If you want to add a title as well please put it in quotes: "Title")
``` 

The authors and links have been added in the description field of the DESCRIPTION file. 

## 2 
``` 
Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. \value{No return value, called for side effects} or similar)
Missing Rd-tags:
     as.matrix.multiPhylo.Rd: \value
     as.vector.phylo.Rd: \value
     plot.troppca.Rd: \value
     read.nexus.to.data.matrix.Rd: \value
     read.tree.to.data.matrix.Rd: \value
     troppca.obj.Rd: \value
``` 

The \\value tags have been added to the to the aforementioned Rd files. 

## 3 

```
\dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ('# Not run:') as a warning for the user.
Does not seem necessary.

Please unwrap the examples if they are executable in < 5 sec, or replace \dontrun{} with \donttest{}.
```

I replaced \\dontrun{} in examples of the following functions to \\donttest{} because they take time > 5 sec:
- `troppca.polytope`
- `troppca.linsp`
- `troppca.linsp2poly` 

## 4 
```
Please make sure that you do not change the user's options, par or working directory. If you really have to do so within functions, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited. e.g.:
...
oldpar <- par(no.readonly = TRUE)    # code line i
on.exit(par(oldpar))            # code line i + 1
...
par(mfrow=c(2,2))            # somewhere after
...
e.g.: R/plot.tropca.R
If you're not familiar with the function, please check ?on.exit. This function makes it possible to restore options before exiting a function even if the function breaks. Therefore it needs to be called immediately after the option change within a function.
``` 
`on.exit()` has been added to `plot.tropca.R` to restore users' settings.


# R CMD check results
There were no ERRORs or WARNINGs. 

There were 3 NOTEs:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Houjie Wang <wanghoujie6688@gmail.com>'

New submission
``` 

This means that this is the first submission of the package to the CRAN. 

```
Possibly misspelled words in DESCRIPTION:
  Yoshida (12:229, 12:349, 12:420)
  Zhang (12:363, 12:433, 12:445)
``` 

These words are the names of some authors of the reference papers describing the methods in the package. 

``` 
Found the following (possibly) invalid DOIs:
  DOI: 10.1093/bioinformatics/btaa564
    From: DESCRIPTION
    Status: Forbidden
    Message: 403
```
It has been checked that this doi is valid.

# Downstream dependencies 
There are currently no downstream dependencies for this package.