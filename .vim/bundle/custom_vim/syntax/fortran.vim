
"colorscheme evening


set foldmethod=syntax
set foldenable " no fold by default
set foldlevel=0   " some value I don't know...
set foldnestmax=10

let fortran_fold=1
let fortran_have_tabs=1
""let fortran_fold_conditionals=1
""let fortran_fold_multilinecomments=1


"hi Mygroup ctermfg=208 ctermbg=black cterm=bold
hi Mygroup ctermfg=208 ctermbg=None cterm=bold
"match Mygroup /^\s*subroutine\w*\|\^\s*function\w*/
"match Mygroup /\s*subroutine\w*/
call matchadd('Mygroup','^\s*subroutine\w*')
call matchadd('Mygroup','^\s*end subroutine\w*')
call matchadd('Mygroup','^\s*function\w*')
call matchadd('Mygroup','^\s*end function\w*')
call matchadd('Mygroup','^\s*module\w*')
call matchadd('Mygroup','^\s*end module\w*')
call matchadd('Mygroup','^\s*program\w*')
call matchadd('Mygroup','^\s*end program\w*')


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


