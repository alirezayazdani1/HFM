help:
	@echo ''
	@echo 'Nektar Utilties:'
	@echo ''
	@echo ' To compile utility run with one of the following suffices:'
	@echo '      dbx ... debugged  version'
	@echo '      opt ... optimised verions'
	@echo ''
	@echo '  i.e. make gridgendbx (generates gridgen with debugging)'
	@echo '       make gridgenopt (generates gridgen optimised)'
	@echo ''
	@echo ' Utilities are:'
	@echo '       gridgen.....Generate mesh points'
	@echo '       nek2tec.....Nektar field to tecplot'
	@echo '       vort........Nektar field to tecplot with vorticity'
	@echo ''
		

%dbx:
	$(MAKE) -f MakeUtils LTYPE=g PROG=$(subst dbx,,$@)

%mopt:
	$(MAKE) -f MakeUtils LTYPE=mopt PROG=$(subst mopt,,$@)

clean:
	rm -f *.o










