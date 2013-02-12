#! bash
#
# Bash completion support for LEMON (commands and --long-options)
#
# Copyright (c) 2013 Victor Terron. All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
# Distributed under the GNU General Public License, version 3.0

# To use these completion routines:
# 1) Copy this file to somewhere (e.g. ~/.lemon-completion.sh).
# 2) Add the following line to your .bashrc/.zshrc:
#     source ~/.lemon-completion.sh


FITS_EXTS="@(fit?(s)|FIT?(S))"

_lemon_import()
{
    local opts
    opts="--object --pattern --counts --filename --follow --exact
    --datek --expk= --objectk --uik"

    if [[ ${cur} != -* ]]; then
        _filedir $FITS_EXTS
    else
	COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
    fi
}

_lemon()
{
    local cur prev commands
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    commands="import"

    # The options that autocomplete depend on the LEMON command being
    # executed. For example, the '--exact' option is specific to the
    # 'import' command, so it must be available only when that is the
    # second word in the current command line (i.e., 'lemon import ...')

    case "${COMP_WORDS[1]}" in
    import)
        _lemon_import
        return 0
	;;
    esac

    COMPREPLY=( $(compgen -W "${commands}" -- ${cur}) )
}

complete -F _lemon lemon
