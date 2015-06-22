echo "Hello Steve!"
export TERM="xterm-color"
export BLASTDB="/usr/local/blastdbs/"
export DYLD_LIBRARY_PATH=/usr/local/mysql/lib/
export HHLIB=/usr/local/bioinf_tools/hhsuite-2.0.16/lib/hh/
export CLASSPATH=.:/usr/local/Java/
export TMPDIR=/Volumes/Zippy/.sysTemp/

#PS1='\[\e[0;33m\]\u\[\e[0m\]@\[\e[0;32m\]\h\[\e[0m\]:\[\e[0;34m\]\w\[\e[0m\]\$ '

alias ls="ls -G"
alias la="ls -a"
alias ll="ls -l"
alias lm="ll | awk '{k=0;for(i=0;i<=8;i++)k+=((substr(\$1,i+2,1)~/[rwx]/)*2^(8-i));if(k)printf(\"%0o \",k);print}'"
export LSCOLORS=gxFxBxDxbxEgEdxbxgxcGx

alias grep="grep --color=auto"

export PATH=$PATH:/usr/local/bioinf:/Users/bondsr/Documents/work_scripts/utilities:/Users/bondsr/Documents/work_scripts/seq_tools:/Users/bondsr/Documents/work_scripts/phylo_tools
export PATH=/usr/local/bin:$PATH

alias ipy="python3 /usr/local/anaconda/bin/ipython notebook"
alias svn=/usr/local/bin/svn-color.py
alias sudo=/usr/local/bin/sudo
source ~/.local/bin/bashmarks.sh


#export PS1="\u@\h "'$(git branch &>/dev/null; if [ $? -eq 0 ]; then \
export PS1="\[\e[0;33m\]\u\[\e[0m\]@\[\e[0;32m\]\h\[\e[0m\]: "'$(git branch &>/dev/null; if [ $? -eq 0 ]; then \
echo "\[\e[0;32m\][\[\e[0;35m\]$(basename `pwd`); \[\e[0;33m\]$(git branch | grep ^*|sed s/\*\ //) \
$(echo `git status` | grep "nothing to commit" > /dev/null 2>&1; if [ "$?" -eq "0" ]; then \
echo -e "\[\e[0;32m\]\xF0\x9F\x8D\xBA  "; else \
echo -e "\[\e[0;31m\]\xF0\x9F\x92\xA9  "; fi)\[\e[0;32m\]]\$"; else \
echo "\[\e[0;35m\][\W]\[\e[m\] \$"; fi) \[\e[0m\]'


# added by Anaconda 1.9.2 installer
export PATH="/usr/local/anaconda/bin:$PATH"

HISTIGNORE="hnote*"
# Used to put notes in a history file
function hnote {
    echo "## NOTE [`date`]: $*" >> $HOME/.history/bash_history-`date +%Y%m`
}

# used to keep my history forever, 'hist' is a link to a python script I wrote (history_recorder.py), that is specific to my own system.
PROMPT_COMMAND="hist \$(history 1) $$"
