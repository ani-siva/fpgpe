default:out

out:
	make -C src run 
	ln -s ./src/run ./out	
clean:
	make -C src clean
full-clean:
	make -C src full-clean
	rm ./out
	rm density/*
