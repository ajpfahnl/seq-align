all: dist

.PHONY: clean dist

clean:
	rm -f *~
	rm -f seq_align.tar.gz
	rm -rf align/__pycache__ align/.DS_Store

setup_files = requirements.txt Makefile setup.sh
run_files = run.sh score_matrix.txt

dist: clean
	tar -zcf seq_align.tar.gz align README.md $(setup_files) $(run_files)



