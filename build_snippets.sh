#!/bin/bash

code_folder=boarding/code
input_code_file=snippets_.c
output_lib_name=snippets

x86_64-apple-darwin13.4.0-clang -I"/Users/andres/opt/miniconda3/lib/R/include" -DNDEBUG -I/Users/andres/opt/miniconda3/lib/R/library/pomp/include -I/Users/andres/research/covid/models/POMP  -D_FORTIFY_SOURCE=2 -mmacosx-version-min=10.9  -I/Users/andres/opt/miniconda3/include  -fPIC  -march=core2 -mtune=haswell -mssse3 -ftree-vectorize -fPIC -fPIE -fstack-protector-strong -O2 -pipe -I/Users/andres/opt/miniconda3/include -fdebug-prefix-map=/opt/concourse/worker/volumes/live/59b7f007-fada-42cf-7435-5bbd0518eaa4/volume/r-base_1570124919999/work=/usr/local/src/conda/r-base-3.6.1 -fdebug-prefix-map=/Users/andres/opt/miniconda3=/usr/local/src/conda-prefix  -c ${code_folder}/${input_code_file} -o ${code_folder}/${output_lib_name}.o

x86_64-apple-darwin13.4.0-clang -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Users/andres/opt/miniconda3/lib/R/lib -Wl,-dead_strip_dylibs -Wl,-pie -Wl,-headerpad_max_install_names -Wl,-dead_strip_dylibs -Wl,-rpath,/Users/andres/opt/miniconda3/lib -L/Users/andres/opt/miniconda3/lib -o ${code_folder}/${output_lib_name}.dylib ${code_folder}/${output_lib_name}.o -L/Users/andres/opt/miniconda3/lib/R/lib -lR -Wl,-framework -Wl,CoreFoundation