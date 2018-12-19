# install the solver 

first, make 2 directories named "build" and "install"  

   mkdir build
   mkdir install  

in the build directory, we have to call the cmake command like this :

   cmake -DUSE_MKL=ON -DCMAKE_INSTALL_PREFIX=<root_dir>/solver/install ..

then, we have to compile the code

   make 

and make the installation 

   make install 

done 
