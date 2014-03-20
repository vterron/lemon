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
XML_EXTS="@(xml|XML)"
LEMONDB_EXTS="@(LEMONdB|lemondb)"

# Match the current word against the list given as argument
_match()
{
    COMPREPLY=( $(compgen -W "${1}" -- ${cur}) )
}

_lemon_import()
{
    local opts
    opts="--object --pattern --counts --filename --follow --exact
    --datek --expk= --objectk --uik"

    if [[ ${cur} != -* ]]; then
        _filedir $FITS_EXTS
    else
	_match "${opts}"
    fi
}

_lemon_seeing()
{
    local opts

    opts="--filename --maximum --margin --snr-percentile --mean
    --sources-percentile --suffix --cores --verbose --fsigma
    --fwhm_dir --esigma --elong_dir --coaddk --saturk --fwhmk"

    if [[ ${cur} != -* ]]; then
        _filedir $FITS_EXTS
    else
	_match "${opts}"
    fi
}

_lemon_mosaic()
{
    local opts checktypes

    opts="--scale --output --overwrite --fraction --name --tempdir
    --check-type --min --max --rak --deck --objectk"

    # The different types of check-image available in SExtractor
    checktypes="NONE IDENTICAL BACKGROUND BACKGROUND_RMS MINIBACK_RMS
    MINIBACKGROUND -BACKGROUND FILTERED OBJECTS -OBJECTS APERTURES
    SEGMENTATION"

    case $prev in
	--output)
	    _filedir $FITS_EXTS
	    return 0
	    ;;
	--check-type)
	    _match "${checktypes}"
	    return 0
	    ;;
    esac

    if [[ ${cur} == -* ]]; then
	_match "${opts}"
    else
        _filedir $XML_EXTS
    fi
}

_lemon_astrometry()
{
    local opts
    opts="--suffix --verbose --ra --dec --radius"
    if [[ ${prev} == --output ]]; then
	_filedir $FITS_EXTS
    elif [[ ${cur} == -* ]]; then
	_match "${opts}"
    else
        _filedir $FITS_EXTS
    fi
}

_lemon_annuli()
{
    local opts
    opts="--output --overwrite --margin --gain --cores --verbose
    --aperture --annulus --dannulus --min-sky --constant
    --minimum-constant --lower --upper --step --sky --width --maximum
    --minimum-images --minimum-stars --pct --weights-threshold
    --max-iters --worst-fraction --expk --coaddk --gaink --uik"

    if [[ ${prev} == --output ]]; then
	_filedir $XML_EXTS
    elif [[ ${cur} == -* ]]; then
	_match "${opts}"
    else
        _filedir $XML_EXTS
    fi
}

_lemon_photometry()
{
    local opts
    opts="--output --overwrite --passband --maximum --margin --gain
    --annuli --cores --verbose --pixels --min-sky --aperture --annulus
    --dannulus --individual-fwhm --aperture-pix --annulus-pix
    --dannulus-pix --expk --coaddk --saturk --gaink --uik"

    case $prev in
	--output)
	    _filedir $LEMONDB_EXTS
	    return 0
	    ;;
	--annuli)
	    _filedir $XML_EXTS
	    return 0
	    ;;
	--pixels)
	    _filedir
	    return 0
	    ;;
    esac

    if [[ ${cur} == -* ]]; then
	_match "${opts}"
    else
        _filedir $XML_EXTS
    fi
}

_lemon_diffphot()
{
    local opts
    opts="--output --overwrite --cores --verbose --minimum-images
    --stars --minimum-stars --pct --weights-threshold --max-iters
    --worst-fraction"

    if [[ ${prev} == --output ]]; then
	_filedir $LEMONDB_EXTS
    elif [[ ${cur} == -* ]]; then
	_match "${opts}"
    else
        _filedir $LEMONDB_EXTS
    fi
}

_lemon_periods()
{

    local opts
    opts="--output --overwrite --initial-step --exhaustive-step
    --cores --verbose"

    if [[ ${prev} == --output ]]; then
	_filedir $LEMONDB_EXTS
    elif [[ ${cur} == -* ]]; then
	_match "${opts}"
    else
        _filedir $LEMONDB_EXTS
    fi
}

_lemon()
{
    local cur prev commands
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    commands="import seeing astrometry mosaic annuli photometry
    diffphot periods juicer"

    # The options that autocomplete depend on the LEMON command being
    # executed. For example, the '--exact' option is specific to the
    # 'import' command, so it must be available only when that is the
    # second word in the current command line (i.e., 'lemon import ...')

    case "${COMP_WORDS[1]}" in
    import)
        _lemon_import
        return 0
	;;
    seeing)
	_lemon_seeing
	return 0
        ;;
    mosaic)
	_lemon_mosaic
	return 0
	;;
    astrometry)
	_lemon_astrometry
	return 0
	;;
    annuli)
	_lemon_annuli
        return 0
	;;
    photometry)
	_lemon_photometry
	return 0
	;;
    diffphot)
	_lemon_diffphot
	return 0
	;;
    periods)
	_lemon_periods
	return 0
	;;
    esac

    _match "${commands}"

}

complete -F _lemon lemon
