#!/bin/bash

# ARDaP Progress Monitor Script
# This script monitors the progress of your ARDaP pipeline

echo "==================================="
echo "ARDaP Pipeline Progress Monitor"
echo "$(date)"
echo "==================================="

# Check if Nextflow is running
echo "1. NEXTFLOW STATUS:"
NEXTFLOW_PID=$(pgrep -f "nextflow.*main.nf")
if [ -n "$NEXTFLOW_PID" ]; then
    echo "   ✓ Nextflow is running (PID: $NEXTFLOW_PID)"
    echo "   Runtime: $(ps -o etime= -p $NEXTFLOW_PID | tr -d ' ')"
else
    echo "   ✗ Nextflow is not running"
fi

echo ""
echo "2. RECENT NEXTFLOW LOG (last 10 lines):"
if [ -f ".nextflow.log" ]; then
    tail -10 .nextflow.log | sed 's/^/   /'
else
    echo "   No .nextflow.log found"
fi

echo ""
echo "3. WORK DIRECTORY ANALYSIS:"
if [ -d "work" ]; then
    TOTAL_TASKS=$(find work -name ".command.run" | wc -l)
    COMPLETED_TASKS=$(find work -name ".exitcode" | wc -l)
    FAILED_TASKS=$(find work -name ".exitcode" -exec cat {} \; | grep -v "^0$" | wc -l)
    RUNNING_TASKS=$((TOTAL_TASKS - COMPLETED_TASKS))
    
    echo "   Total tasks: $TOTAL_TASKS"
    echo "   Completed: $COMPLETED_TASKS"
    echo "   Running: $RUNNING_TASKS"
    echo "   Failed: $FAILED_TASKS"
    
    if [ $FAILED_TASKS -gt 0 ]; then
        echo "   Failed task directories:"
        find work -name ".exitcode" -exec grep -l -v "^0$" {} \; | sed 's/.exitcode//' | sed 's/^/     /'
    fi
else
    echo "   No work directory found"
fi

echo ""
echo "4. CURRENT PROCESSES:"
echo "   ARDaP-related processes:"
ps aux | grep -E "(nextflow|gatk|snpEff|bwa|samtools|trimmomatic)" | grep -v grep | head -10 | awk '{print "     " $2 " " $11 " " $12 " " $13}' || echo "     None found"

echo ""
echo "5. SYSTEM RESOURCES:"
echo "   CPU Usage: $(top -bn1 | grep "Cpu(s)" | awk '{print $2}' | cut -d'%' -f1)%"
echo "   Memory: $(free -m | awk 'NR==2{printf "%.1f%%", $3*100/$2 }')"
echo "   Load Average: $(uptime | awk -F'load average:' '{print $2}')"

echo ""
echo "6. DISK SPACE:"
echo "   Current directory: $(df -h . | awk 'NR==2 {print $4 " free of " $2}')"

echo ""
echo "7. OUTPUTS STATUS:"
if [ -d "Outputs" ]; then
    echo "   Output directories:"
    ls -la Outputs/ | grep "^d" | awk '{print "     " $9 ": " $5 " bytes"}'
else
    echo "   No Outputs directory found yet"
fi

echo ""
echo "==================================="
echo "Monitor completed at $(date)"
echo "Run 'watch -n 30 ./monitor_ardap.sh' for continuous monitoring"
echo "==================================="