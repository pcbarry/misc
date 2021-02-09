set nocompatible              " required
filetype off                  " required

" set the runtime path to include Vundle and initialize
set rtp+=~/.vim/bundle/Vundle.vim
call vundle#begin()

" alternatively, pass a path where Vundle should install plugins
" call vundle#begin('~/some/path/here')

" let Vundle manage Vundle, required
Plugin 'gmarik/Vundle.vim'

" add all your plugins here (note older versions of Vundle
" used Bundle instead of Plugin)

" ...

" All of your Plugins must be added before the following line
call vundle#end()            " required
filetype plugin indent on    " required

" general settings
"-----------------
syntax on 

set backspace=indent,eol,start
set ruler
set number
set ls=2
au BufNewFile,BufRead *.m set filetype=mma
set encoding=utf-8


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
hi Folded ctermfg=100 ctermbg=000
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

" trick for pasting
" -----------------
set pastetoggle=<F2>

" python
" ------
Plugin 'tmhedberg/SimpylFold'
let g:SimpylFold_fold_import = 0
au BufNewFile,BufRead *.py
    \ set tabstop=4 |
    \ set softtabstop=4 |
    \ set shiftwidth=4 |
"    \ set textwidth=79 |
    \ set expandtab |
    \ set autoindent |
    \ set fileformat=unix |

let python_highlight_all=1
Plugin 'vim-scripts/indentpython.vim'
