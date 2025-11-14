edit: ## vim command to initialise editing environment
	vim workflow/Snakefile workflow/rules/* workflow/de-novo-assembly/* 

conon: ##command to activate conda and enviroment
	bash source ~/.profile
	bash conda activate rna-seq
blast_db: ##make blast core_nt database on aws instance
	mkdir blast_db
	aws s3 cp --no-sign-request s3://ncbi-blast-databases/2025-09-11-01-05-02/ blast_db/ --exclude "*" --include "core_nt*" --recursive

#unison setup
PLATFORM := $(shell uname)
ARCH := $(shell uname -m)
UNISON_VERSION := 2.53.7

ifeq ($(PLATFORM),Darwin)
ifeq ($(ARCH),x86_64)
unison_file := "unison-$(UNISON_VERSION)-macos-x86_64.tar.gz"
endif
ifeq ($(ARCH),arm64)
unison_file := "unison-$(UNISON_VERSION)-macos-arm64.tar.gz"
endif
endif
ifeq ($(PLATFORM),Linux)
ifeq ($(ARCH),x86_64)
unison_file := "unison-$(UNISON_VERSION)-ubuntu-x86_64.tar.gz"
endif
ifeq ($(ARCH),arm64)
unison_file := "unison-$(UNISON_VERSION)-ubuntu-arm64.tar.gz"
endif
endif
unison: ##use this to sync repo with a remote host. (this command just installs unison)
	rm -rf bin/unison
	mkdir -p bin/unison
	##download file
	wget --no-check-certificate --content-disposition -P ./bin/unison/ "https://github.com/bcpierce00/unison/releases/download/v$(UNISON_VERSION)/$(unison_file)"
	##extract and tidy
	tar -xvzf bin/unison/$(unison_file) -C bin/unison
	mkdir bin/temp
	mv bin/unison/bin/unison* bin/temp/
	rm -r bin/unison
	mv bin/temp/unison* bin/
	rm -r bin/temp
	#make sure config file is in right place
	#mkdir -p $HOME/.unison
	#cp config/unison/heterogeneous-graphlet-counting.prf $HOME/.unison/
	#set up share directory
	mkdir -p output/share
	# if on mac, set up fsmonitor
	@if [ "$$(uname)" = "Darwin" ]; then \
		brew install autozimu/homebrew-formulas/unison-fsmonitor; \
	fi
env:
	mamba env create --file environment.yml --force
