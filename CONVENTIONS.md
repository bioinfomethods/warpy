# Warpy Development Conventions

## Bash Script Best Practices

### Path Handling

1. **Always use absolute paths for critical directories**
   - When a script needs to reference the project base directory, calculate it relative to the script's location
   - Use `SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"` to get the script's directory
   - Calculate base directory as `WARPY_BASE="$( cd "$SCRIPT_DIR/.." && pwd )"`

2. **Convert relative paths to absolute when accepting user input**
   - If a script accepts file/directory paths via command-line options, convert relative paths to absolute
   - Example:
     ```bash
     if [[ "$OPTARG" = /* ]]; then
       export data_path=$OPTARG
     else
       export data_path="$(cd "$(dirname "$OPTARG")" && pwd)/$(basename "$OPTARG")"
     fi
     ```

3. **Output directories should be relative to project base, not script location**
   - Don't create output directories inside the `scripts/` directory
   - Use the project base directory to determine where outputs should go
   - Example: `BATCH_DIR=$WARPY_BASE/batches/${sample_id}` not `BATCH_DIR=$SCRIPT_DIR/batches/${sample_id}`

### Environment Detection

1. **Support environment variable overrides**
   - Check for environment-specific variables (like `BPIPE_DEFAULT_ENVIRONMENT`) first
   - Fall back to auto-detection only if the variable is not set

2. **Auto-detect platform when possible**
   - Use `$OSTYPE` to detect macOS (`darwin*`)
   - Use `command -v` to check for platform-specific commands (e.g., `sbatch` for HPC environments)
   - Provide sensible defaults for common platforms

3. **Order of precedence for environment detection**
   - Explicit environment variable (highest priority)
   - Platform-specific detection (macOS, HPC with sbatch, etc.)
   - Default fallback (lowest priority)

### Directory Creation

1. **Create parent directories before child directories**
   - Use `mkdir -p` to ensure parent directories exist
   - Example: `mkdir -p "$WARPY_BASE/batches"` before creating `$BATCH_DIR`

### Script Portability

1. **Make scripts work from any directory**
   - Don't assume the script is run from a specific location
   - Calculate all paths relative to the script's actual location
   - Test scripts by running them from different directories
