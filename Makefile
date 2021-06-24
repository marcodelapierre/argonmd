.PHONY: all
all: serial_one serial_multi


.PHONY: serial_one
serial_one:
	cd serial_one_source && $(MAKE)

.PHONY: serial_multi
serial_multi:
	cd serial_multi_sources && $(MAKE)


.PHONY: clean_all
clean_all:
	rm -f */*.x
