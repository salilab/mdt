MODINSTALL=/synth1/home/ben/modeller/current-opt/
PYTHON=python

_mdt.so: mdt_wrap.c
	${PYTHON} setup.py build && cp build/lib*/_mdt.so .

mdt_wrap.c: mdt.i
	swig -python -noproxy $<

clean:
	rm -rf build *.o mdt_wrap.c *.pyc *.so
