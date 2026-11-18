# claude_jobs_guard.ps1 - keeps Claude's MATLAB jobs out of the way of the behaviour experiment (user, 2026-11-17).
# Policy (what the watch does while the task runs; stored in <StateDir>\policy.txt, changed with -Mode policy):
#   freeze  FREEZE (suspend - never kill) Claude's jobs while the task runs. Default.
#   idle    let them run at Idle priority (the experiment programs always get the CPU first). Chosen by the user on
#           2026-11-17 ("네, 지금 시작") after the day-2 sessions that overlapped analysis showed no timing or frame
#           change. A failed per-session check switches back to freeze.
# Modes
#   watch   Run it as its own powershell process, in the background, all day:
#             powershell -NoProfile -ExecutionPolicy Bypass -File <this file> -Mode watch
#           Every few seconds it puts Claude's jobs at Idle priority, forgets stale entries, writes a heartbeat (with
#           the last time the task was seen and the policy) and, under the freeze policy, freezes the jobs while the
#           task (python/pythonw with main_hansol, main.py or main_saleae) runs. It never un-freezes anything.
#   policy  -Policy freeze|idle : set the policy the running watch uses (from its next round).
#   freeze  Freeze Claude's jobs now.
#   resume  Un-freeze Claude's jobs. Needs a watching guard. Under the freeze policy it also refuses while the task
#           runs and within QuietMin minutes of the last session unless -DayOver (= the user said the day is over).
#   ready   Exit 0 only if a guard is watching and no Claude job exists (frozen ones need -Mode resume); under the
#           freeze policy also: the task is not running and the last session ended QuietMin minutes ago (or -DayOver).
#   quiet   Exit 0 only if a guard is watching and, under the freeze policy, the task is not running.
#   status  Show the task, the guard, the policy and Claude's jobs.
# Claude's jobs = matlab.exe / MATLAB.exe started with -batch whose command line has teach_finetune,
# run_all_config_teach, eye_finetune or eye_track (the 9670f black whole-eye model), plus their child processes (never the shared MathWorksServiceHost service). Frozen jobs use
# no CPU and lose no work. The user's own MATLAB desktop (started with -r, or from the Start menu) is never touched.
# History: replaces pause_on_experiment.ps1, which KILLED the jobs; the bash script then started the day-3 job during
# the experiment (2026-11-17 10:49:51-10:54:02). A read-only review (2026-11-17) found and this version fixes: the
# PS 5.1 single-CimInstance .Count bug (the first draft never froze anything), the gap between mice, PID reuse, the
# shared MathWorks service, and a desktop MATLAB started with a teach config.
param(
    [ValidateSet('watch', 'policy', 'freeze', 'resume', 'ready', 'quiet', 'status')][string]$Mode = 'status',
    [ValidateSet('freeze', 'idle')][string]$Policy = 'freeze',
    [switch]$DayOver,
    [switch]$Once,
    [int]$Every = 5,
    [double]$QuietMin = 30,
    [string]$JobWql = "Name='matlab.exe'",
    [string]$JobName = '^matlab\.exe$',
    [string]$JobCmd = '-batch\b.*(teach_finetune|run_all_config_teach|eye_finetune|eye_track)',
    [string]$TaskWql = "Name='python.exe' OR Name='pythonw.exe'",
    [string]$TaskName = '^pythonw?\.exe$',
    [string]$TaskCmd = 'main_hansol|\bmain\.py|main_saleae',
    [string]$StateDir = '',
    [string]$Note = ''
)
if (-not $StateDir) { $StateDir = Join-Path $env:LOCALAPPDATA 'claude_jobs_guard' }   # local, not synced
# note file next to the project's analysis folder (this script lives in <project>\tools); no absolute path stored
if (-not $Note) { $Note = Join-Path (Split-Path $PSScriptRoot -Parent) 'analysis\pupil_teach\models\paused_for_experiment.txt' }
New-Item -ItemType Directory -Force -Path $StateDir | Out-Null
$StateFile = Join-Path $StateDir 'frozen.txt'      # one "pid:creation-ticks" per line
$BeatFile = Join-Path $StateDir 'heartbeat.txt'    # "guard-pid now-ticks task-last-seen-ticks policy"
$PolicyFile = Join-Path $StateDir 'policy.txt'     # "freeze" or "idle"

if (-not ('ClaudeGuard.Nt' -as [type])) {
    Add-Type -Namespace ClaudeGuard -Name Nt -MemberDefinition @'
[System.Runtime.InteropServices.DllImport("ntdll.dll")] public static extern int NtSuspendProcess(System.IntPtr processHandle);
[System.Runtime.InteropServices.DllImport("ntdll.dll")] public static extern int NtResumeProcess(System.IntPtr processHandle);
'@
}

# NOTE: functions unroll arrays, so a single result comes back as a bare CimInstance, which has no .Count in
# Windows PowerShell 5.1. Always call them as @(Get-Task) / @(Get-JobTree) / @(Read-State).
function Get-Task {
    Get-CimInstance Win32_Process -Filter $TaskWql | Where-Object { ($_.Name -match $TaskName) -and ([string]$_.CommandLine -match $TaskCmd) }
}

function Get-JobTree {
    # the matching -batch roots plus their descendants; a child must be younger than its parent (PID reuse), and the
    # shared MathWorksServiceHost service (used by every MATLAB, the user's too) is never included
    $roots = @(Get-CimInstance Win32_Process -Filter $JobWql | Where-Object { ($_.Name -match $JobName) -and ([string]$_.CommandLine -match $JobCmd) })
    if ($roots.Count -eq 0) { return }
    $all = @(Get-CimInstance Win32_Process)
    $byId = @{}
    foreach ($p in $all) { $byId[[int]$p.ProcessId] = $p }
    $ids = New-Object 'System.Collections.Generic.HashSet[int]'
    foreach ($r in $roots) { if ($byId.ContainsKey([int]$r.ProcessId)) { [void]$ids.Add([int]$r.ProcessId) } }
    do {
        $grew = $false
        foreach ($p in $all) {
            $id = [int]$p.ProcessId; $pp = [int]$p.ParentProcessId
            if ($ids.Contains($id) -or (-not $ids.Contains($pp))) { continue }
            if (($null -eq $p.CreationDate) -or ($p.CreationDate -lt $byId[$pp].CreationDate)) { continue }
            if (($p.Name -match '^MathWorksServiceHost') -and ([string]$p.CommandLine -match '\sservice\s|--detached')) { continue }
            [void]$ids.Add($id); $grew = $true
        }
    } while ($grew)
    $all | Where-Object { $ids.Contains([int]$_.ProcessId) }
}

function Get-Key($p) { '{0}:{1}' -f [int]$p.ProcessId, $p.CreationDate.ToUniversalTime().Ticks }

function Test-Frozen([int]$id) {
    # true when every thread of the process waits with reason Suspended
    try {
        $q = [System.Diagnostics.Process]::GetProcessById($id)
        $n = 0
        foreach ($t in $q.Threads) {
            $n++
            if ($t.ThreadState -ne [System.Diagnostics.ThreadState]::Wait) { return $false }
            if ($t.WaitReason -ne [System.Diagnostics.ThreadWaitReason]::Suspended) { return $false }
        }
        return ($n -gt 0)
    } catch { return $false }
}

function Read-State {
    if (Test-Path -LiteralPath $StateFile) { Get-Content -LiteralPath $StateFile | Where-Object { $_ -match '^\d+:\d+$' } }
}

function Write-State($keys) {
    $keys = @($keys | Where-Object { $_ } | Sort-Object -Unique)
    if ($keys.Count -eq 0) { Remove-Item -LiteralPath $StateFile -ErrorAction SilentlyContinue; return }
    $tmp = "$StateFile.tmp"
    Set-Content -LiteralPath $tmp -Value $keys -Encoding Ascii
    [System.IO.File]::Copy($tmp, $StateFile, $true)
    Remove-Item -LiteralPath $tmp -ErrorAction SilentlyContinue
}

function Say([string]$msg) {
    try { Add-Content -LiteralPath $Note -Value $msg -Encoding UTF8 -ErrorAction Stop } catch { }
    Write-Host $msg
}

function Get-StoredPolicy {
    # the policy file wins over the -Policy argument of the watch
    try { $s = ([string](Get-Content -LiteralPath $PolicyFile -ErrorAction Stop | Select-Object -First 1)).Trim() } catch { $s = '' }
    if ($s -eq 'idle' -or $s -eq 'freeze') { return $s }
    return $Policy
}

function Get-Beat {
    # the live guard's heartbeat, or $null: Pid, Age (s since its last round), SinceTask (min since the task was
    # seen), Policy (the one that guard is using)
    try { $b = @(((Get-Content -LiteralPath $BeatFile -ErrorAction Stop) | Select-Object -First 1) -split ' ') } catch { return $null }
    if ($b.Count -lt 3) { return $null }
    $g = @(Get-CimInstance Win32_Process -Filter ('ProcessId=' + [int]$b[0]))
    if (($g.Count -eq 0) -or ([string]$g[0].CommandLine -notmatch 'claude_jobs_guard\.ps1') -or
        ([string]$g[0].CommandLine -notmatch '-Mode\s+watch')) { return $null }
    $now = (Get-Date).Ticks
    $pol = 'freeze'; if (($b.Count -ge 4) -and ($b[3] -eq 'idle')) { $pol = 'idle' }
    [pscustomobject]@{ Pid = [int]$b[0]; Age = ($now - [long]$b[1]) / 1e7; SinceTask = ($now - [long]$b[2]) / 6e8; Policy = $pol }
}

function Test-Live($b) { return (($null -ne $b) -and ($b.Age -le 3 * $Every + 5)) }

function Set-Idle($tree) {
    foreach ($p in $tree) {
        try {
            $q = [System.Diagnostics.Process]::GetProcessById([int]$p.ProcessId)
            if ($q.PriorityClass -ne [System.Diagnostics.ProcessPriorityClass]::Idle) { $q.PriorityClass = 'Idle' }
            $q.Dispose()
        } catch { }
    }
}

function Freeze-Jobs($tree, [string]$why) {
    $state = @(Read-State)
    $new = @(); $adopted = @()
    foreach ($p in $tree) {
        $id = [int]$p.ProcessId; $key = Get-Key $p
        if ($state -contains $key) { continue }
        if (Test-Frozen $id) { $adopted += $key; continue }      # frozen earlier: never suspend twice
        try {
            $q = [System.Diagnostics.Process]::GetProcessById($id)
            if ([ClaudeGuard.Nt]::NtSuspendProcess($q.Handle) -eq 0) { $new += $key } else { Say "freeze FAILED: pid $id" }
        } catch { Say "freeze FAILED: pid $id - $($_.Exception.Message)" }
    }
    if (($new.Count + $adopted.Count) -gt 0) {
        Write-State ($state + $new + $adopted)
        Say ("{0} froze {1} process(es) [{2}] - {3}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), $new.Count,
            ((@($new) + @($adopted) | ForEach-Object { ($_ -split ':')[0] }) -join ' '), $why)
    }
    return $new.Count
}

function Resume-Jobs {
    # every process of the job tree that is listed or frozen (works even if the state file was lost)
    $state = @(Read-State); $done = @()
    foreach ($p in @(Get-JobTree)) {
        $id = [int]$p.ProcessId
        if (-not (($state -contains (Get-Key $p)) -or (Test-Frozen $id))) { continue }
        try {
            $q = [System.Diagnostics.Process]::GetProcessById($id)
            [void][ClaudeGuard.Nt]::NtResumeProcess($q.Handle)
            for ($k = 0; ($k -lt 5) -and (Test-Frozen $id); $k++) { [void][ClaudeGuard.Nt]::NtResumeProcess($q.Handle) }
            $done += $id
        } catch { Say "resume FAILED: pid $id - $($_.Exception.Message)" }
    }
    Write-State @()
    Say ("{0} resumed {1} process(es) [{2}]" -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), $done.Count, ($done -join ' '))
}

function Show-Status {
    $t = @(Get-Task)
    if ($t.Count -gt 0) { Write-Output ('task RUNNING: ' + (($t | ForEach-Object { "$($_.ProcessId) $([string]$_.CommandLine)" }) -join ' | ')) }
    else { Write-Output 'task not running' }
    $b = Get-Beat
    if (-not (Test-Live $b)) { Write-Output ('guard: NOT watching (policy file: {0})' -f (Get-StoredPolicy)) }
    else { Write-Output ('guard: pid {0}, policy {1}, last round {2:N0} s ago, task last seen {3:N1} min ago' -f $b.Pid, $b.Policy, $b.Age, $b.SinceTask) }
    $state = @(Read-State)
    $tree = @(Get-JobTree)
    if ($tree.Count -eq 0) { Write-Output 'no Claude job' }
    foreach ($p in $tree) {
        $id = [int]$p.ProcessId
        $q = $null; try { $q = [System.Diagnostics.Process]::GetProcessById($id) } catch { }
        $cpu = $(if ($q) { $q.TotalProcessorTime.TotalSeconds } else { 0 })
        $pri = $(if ($q) { [string]$q.PriorityClass } else { '?' })
        $st = $(if (Test-Frozen $id) { 'FROZEN' } else { 'running' })
        Write-Output ('{0,6} {1,-22} {2,-7} {3,-11} listed:{4,-5} cpu {5,9:N1} s' -f $id, $p.Name, $st, $pri, ($state -contains (Get-Key $p)), $cpu)
    }
}

switch ($Mode) {
    'status' { Show-Status }
    'policy' {
        Set-Content -LiteralPath $PolicyFile -Value $Policy -Encoding Ascii
        Say ('{0} policy set to {1}' -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), $Policy)
    }
    'quiet' {
        $b = Get-Beat
        if (-not (Test-Live $b)) { Write-Output 'NOT QUIET - no guard is watching'; exit 1 }
        if (($b.Policy -eq 'freeze') -and (@(Get-Task).Count -gt 0)) { Write-Output 'NOT QUIET - the task is running'; exit 1 }
        Write-Output ('quiet (policy {0})' -f $b.Policy); exit 0
    }
    'ready' {
        $b = Get-Beat
        if (-not (Test-Live $b)) { Write-Output 'NOT READY - no guard is watching (start one with -Mode watch)'; exit 1 }
        if ($b.Policy -eq 'freeze') {
            $t = @(Get-Task)
            if ($t.Count -gt 0) { Write-Output "NOT READY - the task is running: $([string]$t[0].CommandLine)"; exit 1 }
        }
        $j = @(Get-JobTree)
        if ($j.Count -gt 0) { Write-Output ('NOT READY - Claude jobs already exist (frozen ones: use -Mode resume): ' + (($j | ForEach-Object { $_.ProcessId }) -join ' ')); exit 1 }
        if (($b.Policy -eq 'freeze') -and (-not $DayOver) -and ($b.SinceTask -lt $QuietMin)) {
            Write-Output ('NOT READY - the task ran {0:N1} min ago (less than {1}); -DayOver only after the user says the sessions are over' -f $b.SinceTask, $QuietMin); exit 1
        }
        Write-State @()
        Write-Output ('ready - policy {0}, the guard (pid {1}) is watching, no Claude job exists' -f $b.Policy, $b.Pid)
        exit 0
    }
    'freeze' { [void](Freeze-Jobs @(Get-JobTree) 'frozen by hand'); Show-Status }
    'resume' {
        $b = Get-Beat
        if (-not (Test-Live $b)) { Write-Output 'NOT RESUMED - no guard is watching (start one with -Mode watch first)'; exit 2 }
        if ($b.Policy -eq 'freeze') {
            $t = @(Get-Task)
            if ($t.Count -gt 0) { Write-Output "NOT RESUMED - the task is running: $([string]$t[0].CommandLine)"; exit 2 }
            if ((-not $DayOver) -and ($b.SinceTask -lt $QuietMin)) {
                Write-Output ('NOT RESUMED - the task ran {0:N1} min ago (less than {1}); -DayOver only after the user says the sessions are over' -f $b.SinceTask, $QuietMin); exit 2
            }
        }
        Resume-Jobs; Show-Status
    }
    'watch' {
        $seen = Get-Date      # its own start counts as a sighting: a session may have just ended
        $was = ''
        Write-Host ('{0} watching every {1} s (state in {2})' -f (Get-Date -Format 'HH:mm:ss'), $Every, $StateDir)
        while ($true) {
            $pol = Get-StoredPolicy
            if ($pol -ne $was) { Write-Host ('{0} policy {1}' -f (Get-Date -Format 'HH:mm:ss'), $pol); $was = $pol }
            $tree = @(Get-JobTree)
            if (@(Get-Task).Count -gt 0) {
                $seen = Get-Date
                if (($pol -eq 'freeze') -and ($tree.Count -gt 0)) { [void](Freeze-Jobs $tree 'the task runs') }
            }
            if ($tree.Count -gt 0) { Set-Idle $tree }
            elseif (Test-Path -LiteralPath $StateFile) { Write-State @() }     # no job left: forget old entries
            try { Set-Content -LiteralPath $BeatFile -Value ('{0} {1} {2} {3}' -f $PID, (Get-Date).Ticks, $seen.Ticks, $pol) -Encoding Ascii -ErrorAction Stop } catch { }
            if ($Once) { break }
            Start-Sleep -Seconds $Every
        }
    }
}
