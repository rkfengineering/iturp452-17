if [ ! -d "build" ]; then
    # Create build directory if it does not exist
    mkdir build
else
    # Clean contents of build directory if it does exist
    rm -r build
    mkdir build
fi