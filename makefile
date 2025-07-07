all:
	python -m build #--no-isolation

local:
	python setup.py build_ext -i

sdist:
	python setup.py sdist

check: sdist
	twine check dist/*

tag:
	## change tag in setup.py
	# git tag -a v0.1.2 -m "version 0.1.2"
	# git push origin v0.1.2

upload: sdist
	twine upload dist/*

clean:
	rm -rf build dist *.egg-info