
# custom modifications
conda config --set changeps1 False
parse_git_branch() { git branch --show-current ; }
export PS1="\e[1;37m(\$(parse_git_branch)) \e[0;33m\w \e[1;37m$\e[0m "
echo "Ready to work!"
echo "Remember to switch to a feature branch before starting the development :)"
