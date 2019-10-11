# Read the pyrhon sources and generate rts files
sphinx-apidoc -feM ../marilib/ -o source/api/
# Build the html doc from rts files
make html
