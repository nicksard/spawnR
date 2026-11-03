package com.nicksard.appsleep;

import android.app.Activity;
import android.app.NotificationManager;
import android.content.ComponentName;
import android.content.Intent;
import android.graphics.Typeface;
import android.net.Uri;
import android.os.Bundle;
import android.provider.Settings;
import android.text.TextUtils;
import android.view.ViewGroup;
import android.widget.Button;
import android.widget.LinearLayout;
import android.widget.ScrollView;
import android.widget.Switch;
import android.widget.TextView;

import java.time.LocalDateTime;
import java.util.List;

public class MainActivity extends Activity {
    private LinearLayout content;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        ScrollView scroll = new ScrollView(this);
        scroll.setFitsSystemWindows(true);
        content = Ui.column(this);
        int pad = Ui.dp(this, 16);
        content.setPadding(pad, pad, pad, Ui.dp(this, 48));
        scroll.addView(content);
        setContentView(scroll);
    }

    @Override
    protected void onResume() {
        super.onResume();
        render();
    }

    private void render() {
        content.removeAllViews();

        boolean guardOn = isAccessibilityServiceOn();
        boolean silencerOn = isNotificationAccessOn();

        // --- Setup -----------------------------------------------------------------------
        LinearLayout guard = Ui.card(this);
        TextView gt = Ui.text(this, guardOn ? "✅  Blocking is on" : "⚠️  Blocking is off", 16);
        gt.setTypeface(Typeface.DEFAULT_BOLD);
        guard.addView(gt);
        if (!guardOn) {
            guard.addView(Ui.secondary(this,
                    "AppSleep needs its Accessibility service turned on to notice when you open a "
                    + "sleeping app.\n\n"
                    + "1. Tap \"Open Accessibility settings\", choose AppSleep, and turn it on.\n"
                    + "2. If the switch is greyed out (\"Restricted setting\"), tap \"Open App info\", "
                    + "then the ⋮ menu at the top right → \"Allow restricted settings\", and try step 1 again."));
            LinearLayout buttons = Ui.row(this);
            buttons.addView(button("Open Accessibility settings",
                    () -> startActivity(new Intent(Settings.ACTION_ACCESSIBILITY_SETTINGS))));
            buttons.addView(button("Open App info", () -> startActivity(
                    new Intent(Settings.ACTION_APPLICATION_DETAILS_SETTINGS,
                            Uri.fromParts("package", getPackageName(), null)))));
            guard.addView(buttons);
        }
        content.addView(guard);

        LinearLayout notif = Ui.card(this);
        TextView nt = Ui.text(this, silencerOn
                ? "✅  Notification silencing is on"
                : "🔔  Notification silencing is off (optional)", 16);
        nt.setTypeface(Typeface.DEFAULT_BOLD);
        notif.addView(nt);
        if (!silencerOn) {
            notif.addView(Ui.secondary(this,
                    "Turn on notification access if you also want sleeping apps' notifications "
                    + "cleared away until they wake up."));
            notif.addView(button("Open notification access",
                    () -> startActivity(new Intent(Settings.ACTION_NOTIFICATION_LISTENER_SETTINGS))));
        }
        content.addView(notif);

        // --- Schedules -------------------------------------------------------------------
        content.addView(Ui.heading(this, "Sleep schedules"));
        List<Schedule> schedules = ScheduleStore.load(this);
        if (schedules.isEmpty()) {
            content.addView(Ui.secondary(this,
                    "No schedules yet. Add one to pick which apps sleep, and when."));
        }
        LocalDateTime now = LocalDateTime.now();
        for (Schedule s : schedules) content.addView(scheduleRow(s, now));

        Button add = button("+  Add schedule",
                () -> startActivity(new Intent(this, EditScheduleActivity.class)));
        LinearLayout.LayoutParams lp = new LinearLayout.LayoutParams(
                ViewGroup.LayoutParams.MATCH_PARENT, ViewGroup.LayoutParams.WRAP_CONTENT);
        lp.topMargin = Ui.dp(this, 8);
        content.addView(add, lp);
    }

    private LinearLayout scheduleRow(Schedule s, LocalDateTime now) {
        LinearLayout card = Ui.card(this);
        LinearLayout top = Ui.row(this);
        LinearLayout texts = Ui.column(this);

        TextView name = Ui.text(this, s.name, 17);
        name.setTypeface(Typeface.DEFAULT_BOLD);
        texts.addView(name);
        texts.addView(Ui.secondary(this, s.timeRangeLabel() + "  ·  " + s.daysLabel()));
        int n = s.packages.size();
        String apps = n + (n == 1 ? " app" : " apps");
        if (s.silenceNotifications) apps += "  ·  notifications silenced";
        texts.addView(Ui.secondary(this, apps));
        if (s.isActiveAt(now)) {
            TextView asleep = Ui.text(this, "🌙  Asleep now", 14);
            asleep.setPadding(0, Ui.dp(this, 4), 0, 0);
            texts.addView(asleep);
        }
        top.addView(texts, Ui.weight1());

        Switch toggle = new Switch(this);
        toggle.setChecked(s.enabled);
        toggle.setOnCheckedChangeListener((b, checked) -> {
            Schedule fresh = ScheduleStore.find(this, s.id);
            if (fresh != null) {
                fresh.enabled = checked;
                ScheduleStore.save(this, fresh);
            }
            render();
        });
        top.addView(toggle);
        card.addView(top);

        card.setOnClickListener(v -> startActivity(new Intent(this, EditScheduleActivity.class)
                .putExtra(EditScheduleActivity.EXTRA_ID, s.id)));
        return card;
    }

    private Button button(String label, Runnable onClick) {
        Button b = new Button(this);
        b.setText(label);
        b.setAllCaps(false);
        b.setOnClickListener(v -> onClick.run());
        return b;
    }

    private boolean isAccessibilityServiceOn() {
        String enabled = Settings.Secure.getString(getContentResolver(),
                Settings.Secure.ENABLED_ACCESSIBILITY_SERVICES);
        if (TextUtils.isEmpty(enabled)) return false;
        ComponentName me = new ComponentName(this, SleepGuardService.class);
        for (String part : enabled.split(":")) {
            ComponentName c = ComponentName.unflattenFromString(part);
            if (me.equals(c)) return true;
        }
        return false;
    }

    private boolean isNotificationAccessOn() {
        NotificationManager nm = getSystemService(NotificationManager.class);
        return nm != null && nm.isNotificationListenerAccessGranted(
                new ComponentName(this, NotificationSilencer.class));
    }
}
