#!/bin/sh
set -e

# Prepend "samtools" if the first argument is not an executable
if ! type -- "$1" &> /dev/null; then
	set -- samtools "$@"
fi

exec "$@"
