.PHONY: test test_remote copy_remote

copy_remote:
	rsync -avz --delete ./ $(ssh):/tmp/magma-unittest

test:
	cd tests && magma -n test_all.m

test_remote: copy_remote
	ssh $(ssh) "cd /tmp/magma-unittest && make test"
