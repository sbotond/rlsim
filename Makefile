VERSION	= 1.4
RD		= rlsim-$(VERSION)_amd64
all: rlsim docs

TARG=rlsim

tools/effest:
	@cd tools; make

rlsim:
	@cd src; make; cp $(TARG) ../

docs:
	@cd doc; make; cd ..

clean:
	rm -f $(TARG)
	@cd src; make clean; cd ..
	@cd doc; make clean; cd ..

release: rlsim docs tools/effest
	@mkdir -p $(RD)/bin
	@cp rlsim $(RD)/bin
	@cp tools/effest $(RD)/bin
	@cp tools/plot_rlsim_report $(RD)/bin
	@cp tools/sel $(RD)/bin
	@cp tools/cov_cmp $(RD)/bin
	@cp tools/pb_plot $(RD)/bin
	@cp doc/rlsim_manual.pdf $(RD)
	@cp README.md $(RD)/
	@cp COPYING $(RD)/
	@tar -cvzf $(RD).tar.gz $(RD) && rm -r $(RD)
	@cp $(RD).tar.gz releases/rlsim-latest_amd64.tar.gz
	@mv $(RD).tar.gz releases/
	@echo Relase $(RD) packaged.
