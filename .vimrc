
" general settings
"-----------------
syntax on 

set backspace=indent,eol,start
set ruler
set number
set ls=2
au BufNewFile,BufRead *.m set filetype=mma


" special trick for markdown files
"----------------------------
let _curfile = expand("%:t")
  if  _curfile =~ ".*\.md"
  filetype plugin indent on
endif


" special trick for makefiles
"----------------------------
let _curfile = expand("%:t")
  if _curfile =~ "Makefile" || _curfile =~ "makefile" || _curfile =~ ".*\.mk"
  set noexpandtab
  set tabstop=2
else
  set expandtab
  set tabstop=2
endif

" trick for folding
"-------------------
nnoremap <space> za

" mouse for xterm
"----------------
set mouse=a
"map <ScrollWheelUp> <C-Y>
"map <ScrollWheelDown> <C-E>
""map <xCSI>[62~ <MouseDown>

if &diff
  colorscheme evening
  set diffopt=filler,context:1000000 " filler is default and inserts empty lines for sync
endif


