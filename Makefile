edit: ## vim command to initialise editing environment
	vim workflow/Snakefile workflow/rules/* workflow/de-novo-assembly/* 

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
