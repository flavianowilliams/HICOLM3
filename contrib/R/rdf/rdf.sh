#!/bin/bash

echo "Please, write de first site"
read var1 && export rdf_a1=$var1

echo "Please, write de second site"
read var2 && export rdf_a2=$var2

Rscript rdf.R

