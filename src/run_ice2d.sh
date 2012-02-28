echo "Cleaning"
make clean
rm ./*.dat
rm -rf RUN
echo "Compiling"
make
echo "Running"
./ice2d