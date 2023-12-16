run: ver cnet pytest 
	rm -f ./test

pytest:
	@PYTHONPATH=../data:${PYTHONPATH} pytest -sv

GCC=g++ -o test -I../../orthogonalstreamlines/intersection/

ver:
	$(GCC) test_vertices.cpp
	./test

cnet:
	$(GCC) test_cablenet.cpp
	./test