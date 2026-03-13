#!/usr/bin/env bash
set -euo pipefail

# ==============================================================
# Slurm interactive session launcher
# ==============================================================
# Uso:
#   bash scripts/slurm_interactive.sh
#
# Opcionalmente puedes sobreescribir parametros por variables de entorno:
#   JOB_TIME=08:00:00 JOB_CPUS=8 JOB_MEM=32G bash scripts/slurm_interactive.sh
#
# O pasar argumentos por CLI:
#   bash scripts/slurm_interactive.sh --time 08:00:00 --cpus 8 --mem 32G --gpus 1 --nodes 1 --ntasks 1
# ==============================================================

# ===================== AJUSTA AQUI EL TAMANO DEL JOB =====================
# 1) Tiempo maximo de la sesion (HH:MM:SS)
JOB_TIME="${JOB_TIME:-04:00:00}"

# 2) Numero de CPU por tarea
JOB_CPUS="${JOB_CPUS:-4}"

# 3) Memoria total solicitada (ej. 8G, 32G, 120G)
JOB_MEM="${JOB_MEM:-16G}"

# 4) GPUs solicitadas (0 si no necesitas GPU)
JOB_GPUS="${JOB_GPUS:-0}"

# 5) Numero de nodos y tareas
JOB_NODES="${JOB_NODES:-1}"
JOB_NTASKS="${JOB_NTASKS:-1}"
# ========================================================================

usage() {
  cat <<'EOF'
Uso:
  bash scripts/slurm_interactive.sh [opciones]

Opciones de tamano del job:
  -t, --time HH:MM:SS    Tiempo maximo de la sesion
  -c, --cpus N           CPUs por tarea
  -m, --mem SIZE         Memoria total (ej. 16G, 64G)
  -g, --gpus N           Numero de GPUs (0 para CPU-only)
  -n, --nodes N          Numero de nodos
  -k, --ntasks N         Numero de tareas

Otras opciones:
  -p, --partition NAME   Particion de Slurm
  -a, --account NAME     Cuenta/proyecto
  -q, --qos NAME         QOS
  -C, --constraint EXPR  Constraint del nodo
  -s, --shell PATH       Shell interactiva (default: /bin/bash)
  -h, --help             Mostrar esta ayuda

Precedencia:
  1) Argumentos CLI
  2) Variables de entorno
  3) Defaults del script
EOF
}

require_value() {
  local opt_name="$1"
  local opt_value="${2:-}"
  if [[ -z "${opt_value}" || "${opt_value}" == -* ]]; then
    echo "Error: la opcion ${opt_name} requiere un valor."
    echo
    usage
    exit 1
  fi
}

# Parametros de cluster/proyecto (ajusta segun tu entorno)
JOB_PARTITION="${JOB_PARTITION:-compute}"
JOB_ACCOUNT="${JOB_ACCOUNT:-}"
JOB_QOS="${JOB_QOS:-}"
JOB_CONSTRAINT="${JOB_CONSTRAINT:-}"

# Shell para la sesion interactiva
INTERACTIVE_SHELL="${INTERACTIVE_SHELL:-/bin/bash}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -t|--time)
      require_value "$1" "${2:-}"
      JOB_TIME="$2"
      shift 2
      ;;
    -c|--cpus)
      require_value "$1" "${2:-}"
      JOB_CPUS="$2"
      shift 2
      ;;
    -m|--mem)
      require_value "$1" "${2:-}"
      JOB_MEM="$2"
      shift 2
      ;;
    -g|--gpus)
      require_value "$1" "${2:-}"
      JOB_GPUS="$2"
      shift 2
      ;;
    -n|--nodes)
      require_value "$1" "${2:-}"
      JOB_NODES="$2"
      shift 2
      ;;
    -k|--ntasks)
      require_value "$1" "${2:-}"
      JOB_NTASKS="$2"
      shift 2
      ;;
    -p|--partition)
      require_value "$1" "${2:-}"
      JOB_PARTITION="$2"
      shift 2
      ;;
    -a|--account)
      require_value "$1" "${2:-}"
      JOB_ACCOUNT="$2"
      shift 2
      ;;
    -q|--qos)
      require_value "$1" "${2:-}"
      JOB_QOS="$2"
      shift 2
      ;;
    -C|--constraint)
      require_value "$1" "${2:-}"
      JOB_CONSTRAINT="$2"
      shift 2
      ;;
    -s|--shell)
      require_value "$1" "${2:-}"
      INTERACTIVE_SHELL="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Error: opcion no reconocida '$1'"
      echo
      usage
      exit 1
      ;;
  esac
done

if ! command -v srun >/dev/null 2>&1; then
  echo "Error: 'srun' no esta disponible en PATH."
  exit 1
fi

SRUN_CMD=(
  srun
  --pty
  --nodes="${JOB_NODES}"
  --ntasks="${JOB_NTASKS}"
  --cpus-per-task="${JOB_CPUS}"
  --mem="${JOB_MEM}"
  --time="${JOB_TIME}"
  --partition="${JOB_PARTITION}"
)

# Si JOB_GPUS > 0, agrega la solicitud de GPU
if [[ "${JOB_GPUS}" -gt 0 ]]; then
  SRUN_CMD+=("--gres=gpu:${JOB_GPUS}")
fi

# Opcionales de cluster
if [[ -n "${JOB_ACCOUNT}" ]]; then
  SRUN_CMD+=("--account=${JOB_ACCOUNT}")
fi

if [[ -n "${JOB_QOS}" ]]; then
  SRUN_CMD+=("--qos=${JOB_QOS}")
fi

if [[ -n "${JOB_CONSTRAINT}" ]]; then
  SRUN_CMD+=("--constraint=${JOB_CONSTRAINT}")
fi

SRUN_CMD+=("${INTERACTIVE_SHELL}")

echo "Lanzando sesion interactiva con:"
printf ' %q' "${SRUN_CMD[@]}"
echo

exec "${SRUN_CMD[@]}"
