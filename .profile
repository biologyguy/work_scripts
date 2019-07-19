echo "Hello Steve!"
export TERM="xterm-color"
export TMPDIR=/Volumes/Zippy/.sysTemp/

alias ls="ls -G"
alias la="ls -a"
alias ll="ls -lh"
alias lm="ll | awk '{k=0;for(i=0;i<=8;i++)k+=((substr(\$1,i+2,1)~/[rwx]/)*2^(8-i));if(k)printf(\"%0o \",k);print}'"
export LSCOLORS=gxFxBxDxbxEgEdxbxgxcGx

alias grep="grep --color=auto"

export PATH=/Users/bondsr/Documents/GitRepos/work_scripts/utilities:/usr/local/bin:$PATH

alias ipy="python3 /usr/local/anaconda3/bin/jupyter notebook"
//alias sudo=/usr/local/bin/sudo
//source ~/.local/bin/bashmarks.sh

export PS1="\[\e[0;33m\]\u\[\e[0m\]@\[\e[0;32m\]\h\[\e[0m\]: "'$(git branch &>/dev/null; if [ $? -eq 0 ]; then \
echo "\[\e[0;32m\][\[\e[0;35m\]$(basename `pwd`); \[\e[0;33m\]$(git branch | grep ^*|sed s/\*\ //) \
$(echo `git status` | grep "nothing to commit" > /dev/null 2>&1; if [ "$?" -eq "0" ]; then \
echo -e "\[\e[0;32m\]\xF0\x9F\x8D\xBA  "; else \
echo -e "\[\e[0;31m\]\xF0\x9F\x92\xA9  "; fi)\[\e[0;32m\]]\$"; else \
echo "\[\e[0;35m\][\W]\[\e[m\] \$"; fi) \[\e[0m\]'


HISTIGNORE="hnote*"
# Used to put notes in a history file
function hnote {
    echo "## NOTE [`date`]: $*" >> $HOME/.history/bash_history-`date +%Y%m`
}

# used to keep my history forever, 'hist' is a link to a python script I wrote (history_recorder.py), that is specific to my own system.
PROMPT_COMMAND="hist \$(history 1) $$"
