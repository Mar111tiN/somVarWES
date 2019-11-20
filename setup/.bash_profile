####################### .bash_profile ################
# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi
################# ENV ##################################

#################### PATHS TO STATIC DATA ###########################

#################### ALIAS ##########################
alias dirs='dirs -v'
alias watch_it='watch -n 60 qstat'
alias refresh='. ~/.bash_profile; . ~/.bashrc'

#################### STYLE ##########################
PS1="\[\033[1;34m\](\u@\h)\[\033[0;34m\] \w \$\[\033[0m\] "

######### PATH ################################################
export PATH=$PATH:$HOME/.local/bin:$HOME/utils/shell
