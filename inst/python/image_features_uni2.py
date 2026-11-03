"""UNI2-h features of polygons, for Giotto's calculateImageFeatures().

load() builds the model once (hf-hub:MahmoodLab/UNI2-h, gated: needs a
Hugging Face token with access). embed(tiles, masks) takes tiles at the
image's own resolution, resizes each to 224 px (bicubic; the mask nearest),
and averages the last layer's patch tokens over the polygon's mask, weighted
by the mask's area in each 14 px patch. A polygon without pixels in its tile
gets the mean of the four centre patches.
"""
import numpy as np
import torch
import timm
from PIL import Image

_MODEL = None
_DEV = None
_MEAN = torch.tensor([0.485, 0.456, 0.406]).view(1, 3, 1, 1)
_STD = torch.tensor([0.229, 0.224, 0.225]).view(1, 3, 1, 1)
TILE, PATCH, GRID = 224, 14, 16


def load(device=None, threads=1):
    """threads: torch CPU threads. 1 by default: inside an R session whose
    data.table already runs OpenMP threads, torch's own OpenMP pool crashes
    the process (segfault) when it starts more than one thread."""
    global _MODEL, _DEV
    if _MODEL is None:
        if threads is not None:
            torch.set_num_threads(int(threads))
        if device is None:
            if torch.cuda.is_available():
                device = "cuda"
            elif torch.backends.mps.is_available():
                device = "mps"
            else:
                device = "cpu"
        kw = dict(img_size=224, patch_size=14, depth=24, num_heads=24,
                  init_values=1e-5, embed_dim=1536, mlp_ratio=2.66667 * 2,
                  num_classes=0, no_embed_class=True,
                  mlp_layer=timm.layers.SwiGLUPacked, act_layer=torch.nn.SiLU,
                  reg_tokens=8, dynamic_img_size=True)
        _MODEL = timm.create_model("hf-hub:MahmoodLab/UNI2-h",
                                   pretrained=True, **kw).eval().to(device)
        _DEV = device
    return _DEV


_PATCH_OF = ((np.arange(TILE)[:, None] // PATCH) * GRID
             + (np.arange(TILE)[None, :] // PATCH)).ravel()
_CENTRE4 = np.zeros(GRID * GRID, np.float32)
_CENTRE4[[7 * GRID + 7, 7 * GRID + 8, 8 * GRID + 7, 8 * GRID + 8]] = 1


def embed(tiles, masks, batch_size=32):
    """tiles (n, h, w, 3) uint8, masks (n, h, w) bool -> (n, 1536) float32."""
    load()
    tiles = np.asarray(tiles, dtype=np.uint8)
    masks = np.asarray(masks, dtype=bool)
    n = tiles.shape[0]
    out = np.zeros((n, 1536), np.float32)
    n_prefix = _MODEL.num_prefix_tokens
    for b in range(0, n, batch_size):
        idx = range(b, min(b + batch_size, n))
        x = np.stack([np.asarray(Image.fromarray(tiles[i]).resize(
            (TILE, TILE), Image.BICUBIC)) for i in idx])
        m = [np.asarray(Image.fromarray(masks[i].astype(np.uint8)).resize(
            (TILE, TILE), Image.NEAREST)).astype(bool).ravel() for i in idx]
        x = torch.as_tensor(x, dtype=torch.float32).permute(0, 3, 1, 2) / 255.0
        x = ((x - _MEAN) / _STD).to(_DEV)
        with torch.inference_mode():
            f = _MODEL.forward_features(x).float().cpu().numpy()
        for k, i in enumerate(idx):
            pt = f[k, n_prefix:]
            w = (np.bincount(_PATCH_OF[m[k]], minlength=GRID * GRID)
                 .astype(np.float32) if m[k].any() else _CENTRE4)
            out[i] = (w @ pt) / w.sum()
    return out
