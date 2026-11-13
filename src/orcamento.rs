// Orcamento de memoria das operacoes fora-da-memoria.
use crate::externo::{Opcoes, OrdemCromossomos};
use std::path::PathBuf;

/// Fracao da RAM livre entregue aos registros. O resto fica para os buffers de
/// leitura e escrita, para o buffer de trabalho do motor e para a folga do
/// alocador, que a conta por registro so estima.
const FRACAO_DA_RAM_LIVRE: f64 = 0.5;

/// Orcamento padrao: metade da RAM livre neste momento, com piso de 256 MB.
pub fn memoria_automatica() -> usize {
    let mut sys = sysinfo::System::new();
    sys.refresh_memory();
    let livre = sys.available_memory() as f64;
    ((livre * FRACAO_DA_RAM_LIVRE) as usize).max(256 << 20)
}

/// Monta as opcoes a partir do que o chamador informou. `None` em
/// `memoria_mb` usa [`memoria_automatica`]; `None` em `dir_temp` usa o
/// diretorio temporario do sistema.
pub fn opcoes(memoria_mb: Option<usize>, dir_temp: Option<&str>, ordem: OrdemCromossomos) -> Opcoes {
    Opcoes {
        memoria_bytes: memoria_mb
            .map(|mb| mb.saturating_mul(1 << 20))
            .unwrap_or_else(memoria_automatica),
        dir_temp: dir_temp.map(PathBuf::from),
        ordem,
    }
}
