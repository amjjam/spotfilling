

build:
	$(MAKE) -C submodules install
	$(MAKE) -C src build

install:

clean:
	$(MAKE) -C src clean
	$(MAKE) -C submodules uninstall

uninstall:
