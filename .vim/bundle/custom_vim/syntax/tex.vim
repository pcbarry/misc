syn spell toplevel
set spell spelllang=en_us
set tw=70

hi SpellBad 	ctermfg=black	ctermbg=LightGrey
"hi SpellRare 	ctermfg=black	ctermbg=LightGrey
"hi SpellCap 	ctermfg=black	ctermbg=LightGrey
"hi SpellLocal 	ctermfg=black	ctermbg=LightGrey

set foldmethod=marker
"set foldenable " no fold by default
"set foldlevel=0   " some value I don't know...
"set foldnestmax=10
set foldmarker=<<<,>>>

set foldtext=CustomFoldText()

fu! CustomFoldText()
	"get first non-blank line
	let fs = v:foldstart
	while getline(fs) =~ '^\s*$' | let fs = nextnonblank(fs + 1)
	endwhile
	if fs > v:foldend
    	let line = getline(v:foldstart)
	else
		let line = substitute(getline(fs), '\t', repeat(' ', &tabstop), 'g')
	endif

	let w = winwidth(0) - &foldcolumn - (&number ? 8 : 0)
	let foldSize = 1 + v:foldend - v:foldstart
	let foldSizeStr = " " . foldSize . " lines "
	let foldLevelStr = repeat("+--", v:foldlevel)
	let lineCount = line("$")
	let foldPercentage = printf("[%.1f", (foldSize*1.0)/lineCount*100) . "%] "
	let expansionString = repeat(".", w - strwidth(
		\ foldSizeStr.line.foldLevelStr.foldPercentage))
	return line . expansionString . foldSizeStr . foldPercentage . foldLevelStr
endf












