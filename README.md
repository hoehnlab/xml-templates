# xml-templates

🚧 Readme is under construction. 🚧

### How to use these templates
These templates are designed for use with Dowser 2.4 or later and BEAST2 2.7.7. 

Download the templates you wish to use (or clone the whole repo), and 
provide the path to these templates to Dowser's getTimeTrees functions.

### How to make your own templates
Make your own templates by changing model specifications in existing templates,
or by creating an example XML in BEAUti, and replacing your alignment and trait 
data with the appropriate variables (`${DATA}`, `${DATE}`, etc). There are some 
required variables, such as DATA, some variables that have special treatment in 
Dowser (DATE, TRAIT, EMP_EQFREQS), and any other variable name you would like
to create can be programmatically specified through Dowser's getTimeTrees.

A more thorough explanation of variables and how they are handled is under 
construction.
