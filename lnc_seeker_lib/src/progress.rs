// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 Arne Kutzner and Pok-Son Kim
use std::sync::Arc;
use std::sync::atomic::{AtomicU32, Ordering};
use pyo3::prelude::*;

pub struct ProgressData {
    pub stage: AtomicU32,
    pub current: AtomicU32,
    pub total: AtomicU32,
}

#[pyclass]
pub struct SessionProgress {
    pub data: Arc<ProgressData>,
}

#[pymethods]
impl SessionProgress {
    #[new]
    fn new() -> Self {
        Self {
            data: Arc::new(ProgressData {
                stage: AtomicU32::new(0),
                current: AtomicU32::new(0),
                total: AtomicU32::new(0),
            }),
        }
    }

    fn get_status(&self) -> (u32, u32, u32) {
        (
            self.data.stage.load(Ordering::SeqCst),
            self.data.current.load(Ordering::SeqCst),
            self.data.total.load(Ordering::SeqCst),
        )
    }
}

pub fn get_session_progress(_id: &str) -> Arc<ProgressData> {
    Arc::new(ProgressData {
        stage: AtomicU32::new(0),
        current: AtomicU32::new(0),
        total: AtomicU32::new(0),
    })
}
