# Create the symlinks to the correct version of the shared object

VER = 1

$(OBJ): lib$(OBJ).so

lib$(OBJ).so: lib$(OBJ).so.$(VER)
	ln -sf lib$(OBJ).so.1 lib$(OBJ).so