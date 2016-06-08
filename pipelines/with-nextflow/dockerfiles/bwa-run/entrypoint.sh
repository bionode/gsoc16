#!/bin/sh
set -e

# Prepend "bwa" if the first argument is not an executable
if ! type -- "$1" &> /dev/null; then
	set -- bwa "$@"
fi

exec "$@"
