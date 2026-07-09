"""Shared UI helpers for Mixing Lab pages.

Provides a consistent "step card" — a bordered container with a numbered
subheader — used to visually group the sequential steps of the protocol and
assessment pages.

Two usage styles are supported:

1. Context-manager (preferred for new code)::

       with step_card("1", "Reaction Kinetics"):
           st.write(...)

2. begin/end pair (for wrapping large existing step bodies without having to
   re-indent them under a ``with`` block)::

       box = begin_step("1", "Reaction Kinetics")
       st.write(...)
       end_step(box)

The begin/end helpers work by entering the container's context manually, so
every ``st.*`` call issued afterwards renders inside the bordered card until
``end_step`` is called. If a step ends the script early via ``st.stop()`` the
container is simply left open for that run — Streamlit resets the render
context on the next run, so no cleanup is required.
"""

from __future__ import annotations

import streamlit as st


def _step_label(number: str, title: str, icon: str = "") -> str:
    label = " · ".join(part for part in (str(number).strip(), title) if part)
    return f"{icon} {label}".strip() if icon else label


def reset_container_stack() -> None:
    """Reset Streamlit's active-container stack to the app root.

    ``begin_step`` enters a container manually (without a ``with`` block) so
    large existing step bodies don't need re-indenting. When a step ends the
    script early via ``st.stop()``, Python never runs the container's
    ``__exit__``, and Streamlit does **not** reset the container stack on a
    ``StopException`` (it only does so for ``st.rerun()``). The half-open
    container then leaks into the *next* run: because the main entry script
    runs first on every rerun, the leak corrupts rendering from the very top
    (page title, guards, buttons) — long before the first ``begin_step`` — and
    makes buttons/widgets require a second click.

    Call this once at the very top of each run (see ``Mixing_Lab.py``), before
    any ``st.*`` rendering, to clear a leak from the previous run. Every step
    card is opened at the app's top level, so resetting to root is always safe.
    Guarded so a Streamlit internals change can never break the app.
    """
    try:
        from streamlit.delta_generator_singletons import (
            context_dg_stack,
            get_default_dg_stack_value,
        )
        context_dg_stack.set(get_default_dg_stack_value())
    except Exception:
        pass


def begin_step(number: str, title: str, icon: str = ""):
    """Open a bordered step card and render its numbered subheader.

    Returns the container so it can be closed later with :func:`end_step`.
    Subsequent ``st.*`` calls render inside the card until it is closed.
    """
    reset_container_stack()
    box = st.container(border=True)
    box.__enter__()
    st.subheader(_step_label(number, title, icon))
    return box


def end_step(box) -> None:
    """Close a step card previously opened with :func:`begin_step`."""
    if box is not None:
        box.__exit__(None, None, None)


class step_card:
    """Context-manager form of :func:`begin_step` / :func:`end_step`."""

    def __init__(self, number: str, title: str, icon: str = ""):
        self._number = number
        self._title = title
        self._icon = icon
        self._box = None

    def __enter__(self):
        self._box = begin_step(self._number, self._title, self._icon)
        return self._box

    def __exit__(self, exc_type, exc_val, exc_tb):
        end_step(self._box)
        return False


# Widget value_types that CANNOT be assigned via st.session_state (buttons,
# download/form-submit buttons, file uploaders, chat inputs). Assigning them
# raises StreamlitValueAssignmentNotAllowedError when the widget is recreated.
_TRIGGER_VALUE_TYPES = frozenset({
    "trigger_value",
    "string_trigger_value",
    "chat_input_value",
    "json_trigger_value",
    "file_uploader_state_value",
})


def unsettable_widget_keys() -> set[str]:
    """Return the session_state keys backed by trigger-type widgets.

    The state-preservation loop in ``Mixing_Lab.py`` re-assigns every
    session_state key to itself so widget values survive page navigation. But
    trigger widgets (buttons etc.) forbid session_state assignment, so touching
    their key triggers a self-heal ``st.rerun()`` that DISCARDS the very click
    that caused the run — making every such button need two clicks.

    This returns those keys so the loop can skip them proactively. Detection
    uses the widget metadata registered on the PREVIOUS run (a trigger widget
    is always clicked on the run *after* it rendered), so the key is skipped
    before it can be reassigned. Guarded so a Streamlit internals change can
    never break the app — on failure it returns an empty set and the reactive
    self-heal still applies.
    """
    keys: set[str] = set()
    try:
        from streamlit.runtime.scriptrunner_utils.script_run_context import (
            get_script_run_ctx,
        )
        ctx = get_script_run_ctx()
        if ctx is None:
            return keys
        state = ctx.session_state._state
        for k in list(st.session_state.keys()):
            try:
                wid = state._get_widget_id(k)
                meta = state._new_widget_state.widget_metadata.get(wid)
                if meta is not None and meta.value_type in _TRIGGER_VALUE_TYPES:
                    keys.add(k)
            except Exception:
                continue
    except Exception:
        pass
    return keys
