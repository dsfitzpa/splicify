"""Verify that multiple tool_use blocks in one response dispatch in parallel,
not sequentially. Each fake tool sleeps 0.1 s; sequential = 0.3 s, parallel
= 0.1 s. Test asserts parallel by checking elapsed wall time is closer to
the single-tool latency than the cumulative one.
"""
import asyncio
import sys
import time
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import agent_v2  # noqa: F401 — triggers path shim
from agent_v2 import main_agent
from agent_v2.main_agent import run_main_agent
from agent_v2.explore import ExploreFinding, run_explore_subagent
from splicify_api.agent.agent_tools import AttachmentRegistry


class _Block:
    def __init__(self, type_, **kwargs):
        self.type = type_
        for k, v in kwargs.items():
            setattr(self, k, v)


class _Resp:
    def __init__(self, content, stop_reason="end_turn"):
        self.content = content
        self.stop_reason = stop_reason


class _SeqMessages:
    def __init__(self, responses):
        self._responses = list(responses)

    def create(self, **kwargs):
        return self._responses.pop(0)


class _SeqClient:
    def __init__(self, responses):
        self.messages = _SeqMessages(responses)


class _SlowDispatch:
    """Each call sleeps 100 ms — parallel dispatch ~100 ms total, sequential ~300 ms."""
    def __init__(self):
        self.calls: list[str] = []

    async def __call__(self, name, args, registry):
        self.calls.append(name)
        await asyncio.sleep(0.1)
        return {"ok": True}


def test_main_agent_dispatches_tools_in_parallel():
    reg = AttachmentRegistry()
    tool_uses = [
        _Block("tool_use", id=f"tu_{i}", name="annotate_attachment",
               input={"attachment_id": f"att_{i}"})
        for i in range(3)
    ]
    client = _SeqClient([
        _Resp(tool_uses, stop_reason="tool_use"),
        _Resp([_Block("text", text="done")], stop_reason="end_turn"),
    ])
    dispatch = _SlowDispatch()

    start = time.perf_counter()
    result = asyncio.run(run_main_agent(
        "x", [], reg, plan_md="",
        client=client, dispatch_fn=dispatch, tools=[],
    ))
    elapsed = time.perf_counter() - start

    # 3 calls × 100 ms = 300 ms sequential. Parallel should be ~100-150 ms.
    assert elapsed < 0.25, f"too slow ({elapsed*1000:.0f} ms) — not parallel"
    assert len(dispatch.calls) == 3
    assert result.n_tool_calls == 3


def test_explore_subagent_dispatches_tools_in_parallel():
    reg = AttachmentRegistry()
    reg.register_product("p", "ACGT" * 100)
    tool_uses = [
        _Block("tool_use", id=f"tu_{i}", name="annotate_attachment",
               input={"attachment_id": f"att_{i}"})
        for i in range(3)
    ]
    client = _SeqClient([
        _Resp(tool_uses, stop_reason="tool_use"),
        _Resp([_Block("text", text='{"summary_md": "ok", "key_facts": {}}')],
              stop_reason="end_turn"),
    ])
    dispatch = _SlowDispatch()

    start = time.perf_counter()
    asyncio.run(run_explore_subagent(
        role="part_scout",
        user_message="x",
        registry=reg,
        tools=[],
        system_prompt="...",
        client=client,
        dispatch_fn=dispatch,
    ))
    elapsed = time.perf_counter() - start

    assert elapsed < 0.25, f"explore loop too slow ({elapsed*1000:.0f} ms) — not parallel"
    assert len(dispatch.calls) == 3
