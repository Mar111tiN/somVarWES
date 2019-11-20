
################## ALIAS ###############################
alias go_bih='sshfs -o follow_symlinks xxxxxxxx@med-transfer1.bihealth.org:/fast/users/xxxxxxx ~/mount'\
' -o volname=BIH_HOME -o allow_other,noapplexattr,noappledouble; cd mount; sublime .'
# shortcut to bash_profile
alias bash_it="nano ~/.bash_profile"
alias dirs="dirs -v"
alias refresh='. ~/.bash_profile; . ~/.bashrc'

##################### ENV VARS ########################
# for charite access
export http_proxy=http://proxy.charite.de:8080
export https_proxy=http://proxy.charite.de:8080
export no_proxy=*.charite.de