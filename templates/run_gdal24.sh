#!/bin/sh

eval "$start"
conda activate $environment
conda info|grep "active environment"

$command
$command2