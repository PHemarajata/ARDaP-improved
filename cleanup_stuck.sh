#!/bin/bash

echo "==================================="
echo "ARDaP Cleanup Script"
echo "$(date)"
echo "==================================="

echo "Killing stuck SnpEff processes from October 17th..."

# Kill old stuck SnpEff processes
ps aux | grep -E "snpEff.*Oct17" | grep -v grep | awk '{print $2}' | xargs -r kill -9

echo "Killed $(ps aux | grep -E "snpEff.*Oct17" | grep -v grep | wc -l) stuck SnpEff processes"

echo ""
echo "Checking for remaining stuck processes..."
ps aux | grep -E "(snpEff|Masked_alignment).*Oct17" | grep -v grep

echo ""
echo "Current Nextflow process:"
ps aux | grep "nextflow.*main.nf" | grep -v grep

echo ""
echo "==================================="
echo "Cleanup completed at $(date)"
echo "==================================="